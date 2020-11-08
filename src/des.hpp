/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef WCS_DES_HPP
#define WCS_DES_HPP

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

//----------------------------
// C++ and C state structures
//----------------------------

/// LP state
struct WCS_LP_State
{
  /// Pointer to an ssa object
  std::unique_ptr<wcs::SSA_NRM> m_ssa_ptr;
  /// Pointer to the reaction network
  std::shared_ptr<wcs::Network> m_net_ptr;
  /// Simulation start time (wall clock)
  double m_t_start;

  WCS_LP_State(std::unique_ptr<wcs::SSA_NRM>&& ssa_ptr,
               std::shared_ptr<wcs::Network>& net_ptr);

  WCS_LP_State(WCS_LP_State&& other) = default;
  WCS_LP_State& operator=(WCS_LP_State&& other) = default;

};

WCS_LP_State::WCS_LP_State(std::unique_ptr<wcs::SSA_NRM>&& ssa_ptr,
                           std::shared_ptr<wcs::Network>& net_ptr)
: m_ssa_ptr(std::move(ssa_ptr)), m_net_ptr(net_ptr), m_t_start(0.0)
{}


/// Structure that contains global state variables
struct WCS_Global_State
{
  Config m_cfg; ///< Configuration state made from command-line arguments
  std::vector<WCS_LP_State> m_LP_states; ///< List of per-LP states
};

/// ROSS-facing structure that does not include any c++ member structure.
struct WCS_State {
  size_t m_lp_idx;
};

/// Event message type
struct WCS_Message {
  unsigned reaction;
};


//-----------------------
// LP type and functions
//-----------------------
void wcs_init(WCS_State *s, tw_lp *lp);
void wcs_prerun(WCS_State *s, tw_lp *lp);
void wcs_event(WCS_State *s, tw_bf *bf, WCS_Message *msg, tw_lp *lp);
void wcs_event_reverse(WCS_State *s, tw_bf *bf, WCS_Message *msg, tw_lp *lp);
void wcs_event_commit(WCS_State *s, tw_bf *bf, WCS_Message *msg, tw_lp *lp);
void wcs_final(WCS_State *s, tw_lp *lp);
tw_peid wcs_map(tw_lpid gid);

tw_lptype wcs_LPs[] = {
  {
    (init_f) wcs_init,
    (pre_run_f) wcs_prerun,
    (event_f) wcs_event,
    (revent_f) wcs_event_reverse,
    (commit_f) wcs_event_commit,
    (final_f) wcs_final,
    (map_f) wcs_map,
    sizeof(WCS_State)
  },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0ul },
};


//-----------------
//  ROSS options
//-----------------
static unsigned int nlp_per_pe = 1u;
static char run_id[1024] = "wcs-ross";

const tw_optdef app_opt[] =
{
  TWOPT_GROUP("WCS ROSS"),
  TWOPT_UINT("nlp", nlp_per_pe, "Number of LPs per processor"),
  TWOPT_CHAR("run", run_id, "User supplied run name"),
  TWOPT_END()
};

#endif // defined(WCS_HAS_ROSS)
#endif // WCS_DES_HPP
