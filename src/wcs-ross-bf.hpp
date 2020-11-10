/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef WCS_ROSS_BF_HPP
#define WCS_ROSS_BF_HPP

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_ROSS)
// ############################################################################
//                        The bit flag definitions
// ############################################################################
// Make sure that those that are used in handling in a same event type must not
// overlap. Bit field id ranges from 0 to 31
#define WCS_BF_FWD       0
#define WCS_BF_SCHED     1

#define WCS_SUB_VAR_(_c_, _i_) _c_ ## _i_
#define WCS_BF_(_bf_,_n_) ((_bf_)->WCS_SUB_VAR_(c, _n_))

#endif // defined(WCS_HAS_ROSS)
#endif // WCS_ROSS_BF_HPP
