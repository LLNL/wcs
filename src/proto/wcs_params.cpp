/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include <string>
#include <iostream>
#include <fstream>
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/message.h>
#include "utils/file.hpp"
#include "proto/utils.hpp"
#include "proto/wcs_params.hpp"

namespace wcs {

#if defined(WCS_HAS_METIS)
static void set_metis_options(const wcs_proto::WCS_Params::Partition_Params& cfg,
                       wcs::Metis_Params& mp, bool verbose = false)
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

  if (verbose) {
    mp.print();
  }
}
#endif // defined(WCS_HAS_METIS)

static void set_SSA_options(const wcs_proto::WCS_Params::Simulation_Params& cfg,
                     wcs::SSA_Params& sp, bool verbose = false)
{
  using sim_params = wcs_proto::WCS_Params::Simulation_Params;

  sp.m_seed = cfg.seed();

  sp.m_max_iter = cfg.max_iter();
  sp.m_max_time = cfg.max_time();

  sp.m_is_iter_set = (sp.m_max_iter > 0u);
  sp.m_is_time_set = (sp.m_max_time > 0.0);

  sp.m_method = static_cast<int>(cfg.method());

  sp.m_tracing = (cfg.trajectory() == sim_params::Tracing);
  sp.m_sampling = (cfg.trajectory() == sim_params::Sampling);

  if (cfg.Sampling_Interval_case() == sim_params::kSamplingIntervalIter)
    sp.m_iter_interval = cfg.sampling_interval_iter();
  else if (cfg.Sampling_Interval_case() == sim_params::kSamplingIntervalTime)
    sp.m_time_interval = cfg.sampling_interval_time();

  sp.m_frag_size = cfg.frag_size();
  sp.m_is_frag_size_set = (cfg.frag_size() != 0u);

  //sp.m_infile = model_file;
  sp.m_outfile = cfg.outfile();
  sp.m_gvizfile = cfg.gvizfile();

  if (!sp.m_is_time_set) {
    sp.m_max_time = wcs::max_sim_time;
  }
  if (!sp.m_is_iter_set && sp.m_is_time_set) {
    sp.m_max_iter = std::numeric_limits<decltype(sp.m_max_iter)>::max();
  }
  if ((sp.m_sampling || sp.m_tracing) && !sp.m_is_frag_size_set) {
    sp.m_frag_size = wcs::default_frag_size;
  }

  if (verbose) {
    sp.print();
  }
}

void read_proto_params(const std::string& filename,
                       wcs::SSA_Params& sp, bool verbose)
{
  wcs_proto::WCS_Params::Simulation_Params wcs_sim_setup;
  wcs::read_prototext(filename, false, wcs_sim_setup);

  if (verbose) {
    std::string str;
    google::protobuf::TextFormat::PrintToString(wcs_sim_setup, &str);
    std::cout << "---- Prototext '" << filename << "' read ----" << std::endl
              << str << std::endl;
  }

  set_SSA_options(wcs_sim_setup, sp, verbose);
}

#if defined(WCS_HAS_METIS)
void read_proto_params(const std::string& filename,
                       wcs::Metis_Params mp, bool verbose)
{
  wcs_proto::WCS_Params::Partition_Params wcs_part_setup;
  wcs::read_prototext(filename, false, wcs_part_setup);

  if (verbose) {
    std::string str;
    google::protobuf::TextFormat::PrintToString(wcs_part_setup, &str);
    std::cout << "---- Prototext '" << filename << "' read ----" << std::endl
              << str << std::endl;
  }

  set_metis_options(wcs_part_setup, mp, verbose);
}

void read_proto_params(const std::string& filename,
                       wcs::SSA_Params& sp, 
                       wcs::Metis_Params& mp,
                       bool verbose)
{
  wcs_proto::WCS_Params wcs_all_setup;
  wcs::read_prototext(filename, false, wcs_all_setup);

  if (verbose) {
    std::string str;
    google::protobuf::TextFormat::PrintToString(wcs_all_setup, &str);
    std::cout << "---- Prototext '" << filename << "' read ----" << std::endl
              << str << std::endl;
  }

  set_SSA_options(wcs_all_setup.sim_setup(), sp, verbose);
  set_metis_options(wcs_all_setup.part_setup(), mp, verbose);
}
#endif // defined(WCS_HAS_METIS)

} // end of namespace wcs
