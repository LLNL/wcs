/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef  __WCS_PROTO_UTILS_HPP__
#define  __WCS_PROTO_UTILS_HPP__

#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if !defined(WCS_HAS_PROTOBUF)
#error WCS requires protocol buffer
#endif

#include <string>
#include <iostream>
#include <fstream>
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/message.h>

namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

void pbuf_log_collector(
       google::protobuf::LogLevel level,
       const char* filename,
       int line,
       const std::string& message);

template<typename T>
bool read_prototext(const std::string& file_name, const bool is_binary,
                    T& wcs_proto_params)
{
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  std::ifstream input(file_name, std::ios::in | std::ios::binary);

  google::protobuf::SetLogHandler(pbuf_log_collector);

  if (!input) {
    std::cerr << file_name << ": File not found!" << std::endl;
    return false;
  }
  if (is_binary) {
    if (!wcs_proto_params.ParseFromIstream(&input)) {
      std::cerr << "Failed to parse WCS_Params in binary-formatted input file: "
                << file_name << std::endl;
      return false;
    }
  } else {
    google::protobuf::io::IstreamInputStream istrm(&input);
    if (!google::protobuf::TextFormat::Parse(&istrm, &wcs_proto_params)) {
      std::cerr << "Failed to parse WCS_Params in text-formatted input file: "
                << file_name << std::endl;
      return false;
    }
  }
  return true;
}

/**@}*/
} // end of namespace wcs
#endif //  __WCS_PROTO_UTILS_HPP__
