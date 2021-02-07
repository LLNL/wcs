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

static void set_metis_options(const wcs_proto::WCS_Params::Partition_Params& cfg,
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

static void set_SSA_options(const wcs_proto::WCS_Params::Simulation_Params& cfg,
                     wcs::SSA_Params& sp)
{
  using sim_params = wcs_proto::WCS_Params::Simulation_Params;

  sp.seed = cfg.seed();

  sp.max_iter = cfg.max_iter();
  sp.max_time = cfg.max_time();

  sp.is_iter_set = (sp.max_iter > 0u);
  sp.is_time_set = (sp.max_time > 0.0);

  sp.method = static_cast<int>(cfg.method());

  sp.tracing = (cfg.trajectory() == sim_params::Tracing);
  sp.sampling = (cfg.trajectory() == sim_params::Sampling);

  if (cfg.Sampling_Interval_case() == sim_params::kSamplingIntervalIter)
    sp.iter_interval = cfg.sampling_interval_iter();
  else if (cfg.Sampling_Interval_case() == sim_params::kSamplingIntervalTime)
    sp.time_interval = cfg.sampling_interval_time();

  sp.frag_size = cfg.frag_size();
  sp.is_frag_size_set = (cfg.frag_size() != 0u);

  //sp.infile = model_file;
  sp.outfile = cfg.outfile();
  sp.gvizfile = cfg.gvizfile();

  if (!sp.is_time_set) {
    sp.max_time = wcs::max_sim_time;
  }
  if (!sp.is_iter_set && sp.is_time_set) {
    sp.max_iter = std::numeric_limits<decltype(sp.max_iter)>::max();
  }
  if ((sp.sampling || sp.tracing) && !sp.is_frag_size_set) {
    sp.frag_size = wcs::default_frag_size;
  }
}

void read_proto_params(const std::string& filename,
                       wcs::SSA_Params& sp)
{
  wcs_proto::WCS_Params::Simulation_Params wcs_sim_setup;
  wcs::read_prototext(filename, false, wcs_sim_setup);
  set_SSA_options(wcs_sim_setup, sp);
}

void read_proto_params(const std::string& filename,
                       wcs::Metis_Params mp)
{
  wcs_proto::WCS_Params::Partition_Params wcs_part_setup;
  wcs::read_prototext(filename, false, wcs_part_setup);
  set_metis_options(wcs_part_setup, mp);
}

void read_proto_params(const std::string& filename,
                       wcs::SSA_Params& sp, 
                       wcs::Metis_Params mp)
{
  wcs_proto::WCS_Params wcs_all_setup;
  wcs::read_prototext(filename, false, wcs_all_setup);

  set_SSA_options(wcs_all_setup.sim_setup(), sp);
  set_metis_options(wcs_all_setup.part_setup(), mp);
}

} // end of namespace wcs
