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
 * This code reads a single prototext input or a set of upto three prototext
 * that conforms to the schema `src/proto/wcs_params.proto`. In case of the
 * former, it includes the input for WCS_Params. In case of the latter, it
 * includes the input for each of Simulation_Params, Partition_Params, and
 * DES_Params.
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
#include "params/wcs_params.hpp"

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

} // end of namespace wcs

int main(int argc, char** argv)
{
  wcs::cmd_line_opts cmd;
  bool ok = cmd.parse_cmd_line(argc, argv);
  if (!ok) return EXIT_FAILURE;
  if (!cmd.m_is_set) return EXIT_SUCCESS;

  cmd.show();

  if (!cmd.m_all_setup.empty()) {
    wcs_proto::WCS_Params wcs_all_setup;
    wcs::read_prototext(cmd.m_all_setup, false, wcs_all_setup);

    std::string str;
    google::protobuf::TextFormat::PrintToString(wcs_all_setup, &str);
    std::cout << str;
  } else {
    if (!cmd.m_sim_setup.empty()) {
      wcs_proto::WCS_Params::Simulation_Params wcs_sim_setup;
      wcs::read_prototext(cmd.m_sim_setup, false, wcs_sim_setup);
      std::string str;
      google::protobuf::TextFormat::PrintToString(wcs_sim_setup, &str);
      std::cout << str;
    }
    if (!cmd.m_part_setup.empty()) {
      wcs_proto::WCS_Params::Partition_Params wcs_part_setup;
      wcs::read_prototext(cmd.m_part_setup, false, wcs_part_setup);
      std::string str;
      google::protobuf::TextFormat::PrintToString(wcs_part_setup, &str);
      std::cout << str;
    }
    if (!cmd.m_des_setup.empty()) {
      wcs_proto::WCS_Params::DES_Params wcs_des_setup;
      wcs::read_prototext(cmd.m_des_setup, false, wcs_des_setup);
      std::string str;
      google::protobuf::TextFormat::PrintToString(wcs_des_setup, &str);
      std::cout << str;
    }
  }

  google::protobuf::ShutdownProtobufLibrary();

  return EXIT_SUCCESS;
}
