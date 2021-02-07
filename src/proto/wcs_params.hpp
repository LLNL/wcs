/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_PROTO_WCS_PARAMS_HPP__
#define  __WCS_PROTO_WCS_PARAMS_HPP__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if !defined(WCS_HAS_PROTOBUF)
#error WCS requires protocol buffer
#endif

#include <string>
#include "proto/wcs_params.pb.h"
#include "params/ssa_params.hpp"
#include "partition/metis_params.hpp"

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

void read_proto_params(const std::string& filename,
                       wcs::SSA_Params& sp);

void read_proto_params(const std::string& filename,
                       wcs::Metis_Params mp);

void read_proto_params(const std::string& filename,
                       wcs::SSA_Params& sp, 
                       wcs::Metis_Params mp);
/**@}*/
} // end of namespace wcs
#endif //  __WCS_PROTO_WCS_PARAMS_HPP__
