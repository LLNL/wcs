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

#if !defined(WCS_HAS_PROTOBUF)
#error WCS requires protocol buffer
#endif

// This code requires METIS
#if defined(WCS_HAS_METIS)
#include <string>
#include <iostream>
#include <vector>
#include <functional>
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include "proto/wcs_params.pb.h"
#include "proto/utils.hpp"
#include "utils/file.hpp"
#include "utils/write_graphviz.hpp"
#include "wcs_params.hpp"
#include "partition/metis_partition.hpp"
#include "partition/partition_info.hpp"
#include "sim_methods/ssa_nrm.hpp"

void set_metis_options(const wcs_proto::WCS_Params::Partition_Params& cfg,
                       wcs::Metis_Params& mp)
{
  mp.m_nparts = cfg.n_parts();
  mobjtype_et objective = (cfg.cut_obj()? METIS_OBJTYPE_CUT : METIS_OBJTYPE_VOL);
  mctype_et coarsening = (cfg.rm_coarse()? METIS_CTYPE_RM : METIS_CTYPE_SHEM);
  mp.set_options(objective, coarsening, cfg.n_iters(), cfg.seed(), cfg.minconn(),
                 cfg.ufactor(), cfg.dbglvl());
  mp.limit_max_vertex_weight(cfg.ub_vwgt());
  mp.set_ratio_of_vertex_weight_to_size(cfg.vratio());
  mp.m_verbose = cfg.verbose();
  mp.m_outfile = cfg.outfile();
}

/// Partition the given reaction network using Metis
bool initial_partition(const wcs::Metis_Params& mp,
                       std::vector<idx_t>& parts,
                       idx_t& objval)
{
  wcs::Metis_Partition partitioner(mp);
  partitioner.prepare();
  //const auto& map_idx2vd = partitioner.get_map_from_idx_to_desc();
  if (mp.m_verbose) {
    partitioner.print_params();
    partitioner.print_metis_graph(std::cout);
    partitioner.print_adjacency(std::cout);
  }

  bool ret = partitioner.run(parts, objval);
  if (!ret) return false;

  const auto& map_idx2desc = partitioner.get_map_from_idx_to_desc();

  for (wcs::partition_id_t i = 0; i < mp.m_nparts; ++i) {
    mp.m_rnet->set_partition(map_idx2desc, parts, i);

    if (!(mp.m_rnet->my_reaction_list()).empty()) {
      const auto gpart_name
        = wcs::append_to_stem(mp.m_outfile, "-" + std::to_string(i));

      wcs::SSA_NRM nrm(mp.m_rnet);
      nrm.init(1, 0.0, mp.get_seed());
      const auto r = nrm.choose_reaction();
      const auto rname = ((mp.m_rnet->graph())[r.second]).get_label();
      std::cout << "First reaction from partition " << i << " is "
                << rname << " at time "  << r.first << std::endl;

      if (!wcs::write_graphviz(gpart_name, mp.m_rnet->graph(), i)) {
        std::cerr << "Failed to write " << gpart_name << std::endl;
        continue;
      }
      std::cout << "(" + std::to_string((mp.m_rnet->my_reaction_list()).size())
        + " reactions, " + std::to_string((mp.m_rnet->my_species_list()).size())
        + " species)" << std::endl;
    }
  }
  return true;
}


int main(int argc, char** argv)
{
  wcs::cmd_line_opts cmd;
  bool ok = cmd.parse_cmd_line(argc, argv);
  if (!ok) return EXIT_FAILURE;
  if (!cmd.m_is_set) return EXIT_SUCCESS;

  cmd.show();

  wcs::Metis_Params mp;

  if (!cmd.m_all_setup.empty()) {
    wcs_proto::WCS_Params wcs_all_setup;
    wcs::read_prototext(cmd.m_all_setup, false, wcs_all_setup);
    set_metis_options(wcs_all_setup.part_setup(), mp);
  } else {
    if (!cmd.m_sim_setup.empty()) {
      wcs_proto::WCS_Params::Simulation_Params wcs_sim_setup;
      wcs::read_prototext(cmd.m_sim_setup, false, wcs_sim_setup);
    }
    if (!cmd.m_part_setup.empty()) {
      wcs_proto::WCS_Params::Partition_Params wcs_part_setup;
      wcs::read_prototext(cmd.m_part_setup, false, wcs_part_setup);
      set_metis_options(wcs_part_setup, mp);

    }
    if (!cmd.m_des_setup.empty()) {
      wcs_proto::WCS_Params::DES_Params wcs_des_setup;
      wcs::read_prototext(cmd.m_des_setup, false, wcs_des_setup);
    }
  }
  google::protobuf::ShutdownProtobufLibrary();

  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;

  rnet.load(cmd.m_input_model);
  rnet.init();

  mp.m_rnet = rnet_ptr;
  mp.m_infile = cmd.m_input_model;

  std::vector<idx_t> parts; ///< Partition assignment result
  idx_t objval; /// Total comm volume or edge-cut of the solution
  ok = initial_partition(mp, parts, objval);
  if (!ok) return EXIT_FAILURE;

  wcs::Partition_Info pinfo(rnet_ptr);
  pinfo.scan(mp.m_verbose);
  pinfo.report();

  return EXIT_SUCCESS;
}
#else  // defined(WCS_HAS_METIS)
#error This code requires METIS
#endif // defined(WCS_HAS_METIS)
