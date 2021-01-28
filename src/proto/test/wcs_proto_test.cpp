/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

/*
 * This code reads a single prototext input that conforms to the schema
 * `src/proto/wcs_params.proto`. The input file can be either text or
 * binary. This compiles independently of WCS.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <tuple>
#include <list>
#include <functional>
#include <google/protobuf/text_format.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include "proto/wcs_params.pb.h"
#include "wcs_params.hpp"


namespace wcs {

void pbuf_log_collector(
       google::protobuf::LogLevel level,
       const char* filename,
       int line,
       const std::string& message)
{
  std::string errmsg
    = std::to_string(static_cast<int>(level)) + ' ' + std::string{filename}
    + ' ' + std::to_string(line) + ' ' + message;
  std::cerr << errmsg << std::endl;
}

bool read_prototext(const std::string& file_name, const bool is_binary,
                    wcs_proto::WCS_Params& wcs_proto_params)
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

} // end of namespace wcs

int main(int argc, char** argv)
{
  if ((argc < 2) || (argc > 3)) {
    std::cout << "Usage: " << argv[0] << " prototext_file is_binary[0|1]" << std::endl;
    return EXIT_SUCCESS;
  }
  std::string prototext = argv[1];
  const bool is_binary = ((argc == 3)? static_cast<bool>(argv[2]) : false);
  wcs_proto::WCS_Params wcs_proto_params;

  wcs::read_prototext(prototext, is_binary, wcs_proto_params);

  std::string str;
  google::protobuf::TextFormat::PrintToString(wcs_proto_params, &str);
  std::cout << str;

  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}
