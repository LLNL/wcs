/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#include <iostream>
#include <getopt.h>
#include <boost/filesystem.hpp>
#include "utils/write_graphviz.hpp"
#include "partition/partition.hpp"
#include "partition/partition_info.hpp"
#include "sim_methods/ssa_nrm.hpp"
#include "utils/file.hpp"
#include <fstream>

#if defined(WCS_HAS_METIS)
#include "partition/metis_partition.hpp"
#else
using idx_t = int;
#endif // defined(WCS_HAS_METIS)


#define OPTIONS "hb:cd:ei:mp:o:rs:u:v:"
static const struct option longopts[] = {
    {"help",      no_argument,        0, 'h'},
    {"ub_vwgt",   required_argument,  0, 'b'},
    {"cut_obj",   no_argument,        0, 'c'},
    {"dbglvl",    required_argument,  0, 'd'},
    {"embedded",  no_argument,        0, 'e'},
    {"n_iters",   required_argument,  0, 'i'},
    {"minconn",   no_argument,        0, 'm'},
    {"n_parts",   required_argument,  0, 'p'},
    {"outfile",   required_argument,  0, 'o'},
    {"rm_coarse", no_argument,        0, 'r'},
    {"seed",      required_argument,  0, 's'},
    {"ufactor",   required_argument,  0, 'u'},
    {"vratio",    required_argument,  0, 'v'},
    { 0, 0, 0, 0 },
};


struct Config {
  idx_t n_iters; ///< Number of refinement iterations
  idx_t n_parts; ///< Number of partitions
  idx_t seed;    ///< Random number seed
  idx_t ufactor; ///< Uniform load balance factor
  idx_t ub_vwgt; ///< Upper-bound of vertex weight
  double vratio; ///< Ratio of vertex weight to size
  bool rm_coarse; ///< Coarsening by random matching
  bool minconn; ///< Minimize the maximum connectivity
  bool cut_obj; ///< optimize edge-cut rather than communication volume
  bool run_embedded; ///< Whether to run the hard-coded example
  bool verbose; ///< Whether to show details
  idx_t dbglvl; ///< Metis debug level

  std::string infile;
  std::string outfile;

  Config();
  void getopt(int& argc, char** &argv);
  void print() const;
  void print_usage(const std::string exec, int code);
};

Config::Config()
: n_iters(10), n_parts(2), seed(7177), ufactor(300), ub_vwgt(1), vratio(1.0),
  rm_coarse(false), minconn(false), cut_obj(false), run_embedded(false),
  verbose(false), dbglvl(0), infile(""), outfile("")
{}

void Config::getopt(int& argc, char** &argv)
{
  int c = 0;
  while ((c = getopt_long(argc, argv, OPTIONS, longopts, NULL)) != -1) {
    switch (c) {
      case 'h': /* --help */
        print_usage(argv[0], 0);
        break;
      case 'b': /* --ub_vwgt */
        ub_vwgt = static_cast<idx_t>(atoi(optarg));
        break;
      case 'c': /* --cut_obj */
        cut_obj = true;
        break;
      case 'd': /* --dbglvl */
        verbose = (atoi(optarg) >= 512);
        dbglvl = static_cast<idx_t>(atoi(optarg) & 511);
        break;
      case 'e': /* --embedded */
        run_embedded = true;
        break;
      case 'i': /* --n_iters */
        n_iters = static_cast<idx_t>(atoi(optarg));
        if (n_iters < static_cast<idx_t>(1)) {
          std::cerr << "Invalid number of refinement iterations!" << std::endl;
          exit(1);
        }
        break;
      case 'm': /* --minconn */
        minconn = true;
        break;
      case 'p': /* --n_parts */
        n_parts = static_cast<idx_t>(atoi(optarg));
        break;
      case 'o': /* --outfile */
        outfile = std::string(optarg);
        break;
      case 'r': /* --rm_coarse */
        rm_coarse = true;
        break;
      case 's': /* --seed */
        seed = static_cast<idx_t>(atoi(optarg));
        break;
      case 'u': /* --ufactor */
        ufactor = static_cast<idx_t>(atoi(optarg));
        break;
      case 'v': /* --vratio */
        vratio = atof(optarg);
        break;
      default:
        std::cerr << "Invalid option detected!" << std::endl;
        print_usage(argv[0], 1);
        break;
    }
  }

  if (optind != (argc - 1)) {
    print_usage (argv[0], 1);
  }

  infile = argv[optind];

  if (outfile.empty()) {
    boost::filesystem::path path = infile;
    std::string base = path.stem().string();
    outfile = base + ".dot";
  }
}

void Config::print() const
{
  using namespace std;
  cout << endl << "========== Config =============" << endl;
  cout << " - n_iters: " << n_iters << endl;
  cout << " - n_parts: " << n_parts << endl;
  cout << " - seed: " <<  seed << endl;
  cout << " - ufactor: " << ufactor << endl;
  cout << " - ub_vwgt: " <<  ub_vwgt << endl;
  cout << " - vratio: " << vratio << endl;
  cout << " - rm_coarse: " <<  (rm_coarse? "true" : "false") << endl;
  cout << " - minconn: " << (minconn? "true" : "false") << endl;
  cout << " - cut_obj: " << (cut_obj? "true" : "false") << endl;
  cout << " - run_embedded: " << (run_embedded? "true" : "false") << endl;
  cout << " - verbose: " << (verbose? "true" : "false") << endl;
  cout << " - dbglvl: " << dbglvl << endl;
  cout << " - infile: " << infile << endl;
  cout << " - outfile: " << outfile << endl;
  cout << "===============================" << endl << endl;
}

void Config::print_usage(const std::string exec, int code)
{
  std::cerr <<
    "Usage: " << exec << " <filename>.graphml\n"
    "    -h, --help\n"
    "            Display this usage information\n"
    "\n"
    "    -b, --ub_vwgt\n"
    "            Upper-bound on the vertex weight\n"
    "\n"
    "    -c, --cut_obj\n"
    "            Objective is to minimize edge-cut.\n"
    "            By default, minimize communication volume.\n"
    "\n"
    "    -d, --dbglvl\n"
    "            If greater than or equal to 512, verbosity is set\n"
    "            The 9 least significant bits are used for Metis debug level.\n"
    "\n"
    "    -e, --embedded\n"
    "            Run a hard-coded example without using an input graph\n"
    "            The input graph file is ignored.\n"
    "\n"
    "    -i, --n_iters\n"
    "            Number of refinement iterations. 10 by default\n"
    "\n"
    "    -m, --minconn\n"
    "            Minimize the maximum connectivity\n"
    "\n"
    "    -p, --n_parts\n"
    "            Number of partitions\n"
    "\n"
    "    -o, --outfile\n"
    "            Specify the output file name\n"
    "\n"
    "    -r, --rm_coarse\n"
    "            Coarsening by random matching.\n"
    "            Sorted heavy-edge matching by default.\n"
    "\n"
    "    -s, --seed\n"
    "            Random number generator seed\n"
    "\n"
    "    -u, --ufactor\n"
    "            Maximum allowed load imbalance factor.\n"
    "\n"
    "    -v, --vratio\n"
    "            Ratio of the vertex weight to the vertex size. (default 1.0)\n"
    "            It is used to generate vertex sizes based on vertex weights.\n"
    "\n";
  exit(code);
}


#if defined(WCS_HAS_METIS)
//-----------------------------------------------------------------------------
void set_metis_options(const Config& cfg, wcs::Metis_Params& mp)
{
  mobjtype_et objective = (cfg.cut_obj? METIS_OBJTYPE_CUT : METIS_OBJTYPE_VOL);
  mctype_et coarsening = (cfg.rm_coarse? METIS_CTYPE_RM : METIS_CTYPE_SHEM);
  mp.set_options(objective, coarsening, cfg.n_iters, cfg.seed, cfg.minconn,
                 cfg.ufactor, cfg.dbglvl);
  mp.limit_max_vertex_weight(cfg.ub_vwgt);
  mp.set_ratio_of_vertex_weight_to_size(cfg.vratio);
}


/// Run Metis on a manually constructed adjacency list
bool run_metis(wcs::Metis_Params& mp, std::vector<idx_t>& parts,
              idx_t& objval, bool verbose = false)
{
  if (verbose) {
    mp.print();
  }

  // Indexes of starting points in adjacent array
  //std::vector<idx_t> xadj = {0, 2, 5, 7, 9, 12, 14, 15};
  std::vector<idx_t> xadj = {0, 2, 5, 7, 9, 12, 15, 16};

  // Adjacent vertices in consecutive index order
  //std::vector<idx_t> adjncy = {1, 3, 0, 2, 4, 1, 5, 0, 4, 1, 3, 5, 2, 6, 5};
  std::vector<idx_t> adjncy   = {1, 3, 0, 2, 4, 1, 5, 0, 4, 1, 3, 5, 2, 4, 6, 5};

  idx_t nvtxs = xadj.size() - 1;
  idx_t nvwghts = 1;

  idx_t max_vwgt = (mp.m_vwgt_max > 0)? mp.m_vwgt_max : 20;

  if (max_vwgt/2 == static_cast<idx_t>(0)) {
    std::cerr << "This problem requires a maximum vertex weight larger than 1!"
              << std::endl;
    return false;
  }

  // vertex weights
  //std::vector<idx_t> vwgt = {1, 10, 1, 10, 1, 20, 1};
  std::vector<idx_t> vwgt = {1, max_vwgt/2, 1, max_vwgt/2, 1, max_vwgt, 1};

  // vertex sizes
  std::vector<idx_t> vsize(vwgt);
  if ((mp.m_ratio_w2s > 0.0) && (mp.m_ratio_w2s != 1.0)) {
    const double ratio_w2s = mp.m_ratio_w2s;
    std::for_each(vsize.begin(), vsize.end(),
      [ratio_w2s](idx_t& s) {
        s = (s <= ratio_w2s)?
              static_cast<idx_t>(1) :
              static_cast<idx_t>(s/ratio_w2s);
      });
  }

  wcs::Metis_Partition::print_metis_inputs(xadj, adjncy, vwgt, vsize, std::cout);

  parts.clear();
  parts.resize(nvtxs);

  int ret
    = METIS_PartGraphKway(&nvtxs, &nvwghts, xadj.data(), adjncy.data(),
                          vwgt.data(), vsize.data(), NULL, &(mp.m_nparts), NULL,
                          NULL, mp.m_opts.data(), &objval, parts.data());

  return wcs::Metis_Partition::check_run(ret, true);
}


/// Partition the given reaction network using Metis
bool initial_partition(const Config& cfg,
                       const std::shared_ptr<wcs::Network>& rnet_ptr,
                       std::vector<idx_t>& parts,
                       idx_t& objval)
{
  wcs::Metis_Params mp;

  mp.set(cfg.n_parts, rnet_ptr);
  set_metis_options(cfg, mp);

  wcs::Metis_Partition partitioner(mp);
  partitioner.prepare();
  //const auto& map_idx2vd = partitioner.get_map_from_idx_to_desc();
  if (cfg.verbose) {
    partitioner.print_params();
    partitioner.print_metis_graph(std::cout);
    partitioner.print_adjacency(std::cout);
  }

  bool ret = partitioner.run(parts, objval);
  if (!ret) return false;

  const auto& map_idx2desc = partitioner.get_map_from_idx_to_desc();

  for (wcs::partition_id_t i = 0; i < mp.m_nparts; ++i) {
    rnet_ptr->set_partition(map_idx2desc, parts, i);

    if (!(rnet_ptr->my_reaction_list()).empty()) {
      const auto gpart_name
        = wcs::append_to_stem(cfg.outfile, "-" + std::to_string(i));

      wcs::SSA_NRM nrm(rnet_ptr);
      nrm.init(1, 0.0, cfg.seed);
      const auto r = nrm.choose_reaction();
      const auto rname = ((rnet_ptr->graph())[r.second]).get_label();
      std::cout << "First reaction from partition " << i << " is "
                << rname << " at time "  << r.first << std::endl;

      if (!wcs::write_graphviz(gpart_name, rnet_ptr->graph(), i)) {
        std::cerr << "Failed to write " << gpart_name << std::endl;
        continue;
      }
      std::cout << "(" + std::to_string((rnet_ptr->my_reaction_list()).size())
        + " reactions, " + std::to_string((rnet_ptr->my_species_list()).size())
        + " species)" << std::endl;
    }
  }
  return true;
}
//-----------------------------------------------------------------------------
#endif // defined(WCS_HAS_METIS)


int main(int argc, char** argv)
{
  Config cfg;
  cfg.getopt(argc, argv);
  if (cfg.verbose) {
    cfg.print();
  }

  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;
  rnet.load(cfg.infile);
  rnet.init();

  if (!wcs::write_graphviz(cfg.outfile, rnet.graph())) {
    std::cerr << "Failed to write " << cfg.outfile << std::endl;
    return EXIT_FAILURE;
  }

  wcs::Partition pg;
  pg.make_intermediate_graph(rnet);

  const auto intm_gfile_name = wcs::append_to_stem(cfg.outfile, "-intm");
  if (!wcs::write_graphviz(intm_gfile_name, pg.intm_graph())) {
    std::cerr << "Failed to write " << intm_gfile_name << std::endl;
    return EXIT_FAILURE;
  }

#if defined(WCS_HAS_METIS)
  std::vector<idx_t> parts; ///< Partition assignment result
  idx_t objval; /// Total comm volume or edge-cut of the solution

  bool ok = true;

  if (cfg.run_embedded) { // A simple test case
    wcs::Metis_Params mp;
    mp.m_nparts = cfg.n_parts;
    set_metis_options(cfg, mp);
    ok = run_metis(mp, parts, objval, cfg.verbose);
  } else {
    ok = initial_partition(cfg, rnet_ptr, parts, objval);
  }

  if (!ok) {
    std::cerr << "Failed to run Metis" << std::endl;
    return EXIT_FAILURE;
  }

  wcs::Partition_Info pinfo(rnet_ptr);
  pinfo.scan(cfg.verbose);
  pinfo.report(cfg.dbglvl < 513);

  std::cout << "objval = " << objval << std::endl;

  std::vector<size_t> count_part(cfg.n_parts);
  for (size_t i = 0ul; i < parts.size(); ++ i) {
    count_part.at(parts[i]) ++;
  }

  for (size_t i = 0ul; i < count_part.size(); ++i) {
    std::cout << "p" << i << " " << count_part[i] << std::endl;
  }
#endif // defined(WCS_HAS_METIS)

  return EXIT_SUCCESS;
}

