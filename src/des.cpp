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


#if defined(WCS_HAS_ROSS)
#include "sim_methods/ssa_nrm.hpp"
#include "ross.h"
#include <cstdio>
#include <cstring>


struct wcs_state
{
   wcs::SSA_NRM* ssa;
   std::shared_ptr<wcs::Network> net_ptr;
};

typedef enum {
  WCS_INIT,
  TENTATIVE_REACTION,
} wcs_event_type;

typedef struct {
   wcs::Network::v_desc_t fired_reaction;
} wcs_message;

void wcs_init(wcs_state *s, tw_lp *lp);
void wcs_event(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp);
void wcs_event_reverse(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp);
void wcs_event_commit(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp);
void wcs_final(wcs_state *s, tw_lp *lp);

tw_lptype wcs_LPs[] = {
  {
    (init_f) wcs_init,
    (pre_run_f) NULL,
    (event_f) wcs_event,
    (revent_f) wcs_event_reverse,
    (commit_f) wcs_event_commit,
    (final_f) wcs_final,
    (map_f) NULL,
    sizeof(wcs_state)
  },
  { NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0ul },
};

const tw_optdef wcs_opts[] = {
  TWOPT_GROUP("ROSS Model"),
  TWOPT_END(),
};



int main (int argc, char* argv[])
{
  tw_opt_add(wcs_opts);
  tw_init(&argc, &argv);

  int num_LPs_per_pe = 1;

  printf("tw_nnodes %d\n", tw_nnodes());
  printf("num_LPs_per_pe %d\n", num_LPs_per_pe);
  printf("message size %lu\n", sizeof(wcs_message));

  tw_define_lps(num_LPs_per_pe, sizeof(wcs_state));

  std::shared_ptr<wcs::Network> rnet_ptr = std::make_shared<wcs::Network>();
  wcs::Network& rnet = *rnet_ptr;
  rnet.load(argv[1]);
  rnet.init();
  rnet.graph();
  
  g_tw_lp_types = wcs_LPs;
  tw_lp_setup_types();
  tw_end();

  return 0;
}

void wcs_init(wcs_state *s, tw_lp *lp)
{
   s->ssa->init_des(s->net_ptr, lp->gid);
}

void wcs_event(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp)
{
   wcs::SSA_NRM::priority_t firing;
   firing.first = tw_now(lp);
   firing.second = msg->fired_reaction;
   s->ssa->forward_des(firing);

   const auto new_firing = s->ssa->choose_reaction();

   tw_event* next_event = tw_event_new(lp->gid, new_firing.first, lp);
   wcs_message* next_msg = reinterpret_cast<wcs_message*>(tw_event_data(next_event));
   next_msg->fired_reaction = new_firing.second;
   tw_event_send(next_event);
}

void wcs_event_reverse(wcs_state *s, tw_bf *bf, wcs_message *msg, tw_lp *lp)
{
   wcs::SSA_NRM::priority_t firing;
   firing.first = tw_now(lp);
   firing.second = msg->fired_reaction;
   s->ssa->backward_des(firing);
}

void wcs_event_commit(wcs_state* s, tw_bf *bf, wcs_message *msg, tw_lp *lp)
{
   s->ssa->commit_des();
}

void wcs_final(wcs_state *s, tw_lp *lp) 
{
   
}
#endif // defined(WCS_HAS_ROSS)
