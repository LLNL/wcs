/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef WCS_HYBRID_HPP
#define WCS_HYBRID_HPP

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_ROSS)

#include <vector>
#include <memory>
#include <ross.h>
#include "reaction_network/network.hpp"
#include "sim_methods/ssa_nrm.hpp"
#include "params/ssa_params.hpp"



void do_nothing(void *,void *,void *);
tw_peid ctr_map(tw_lpid gid);
void post_lpinit_setup(tw_pe* pe);




// //-----------------
// //  ROSS options
// //-----------------
// static unsigned int nlp_per_pe = 1u;
// static char run_id[1024] = "wcs-ross";

// const tw_optdef app_opt[] =
// {
//   TWOPT_GROUP("WCS ROSS"),
//   TWOPT_UINT("nlp", nlp_per_pe, "Number of LPs per processor"),
//   TWOPT_CHAR("run", run_id, "User supplied run name"),
//   TWOPT_END()
// };

#endif // defined(WCS_HAS_ROSS)
#endif // WCS_HYBRID_HPP
