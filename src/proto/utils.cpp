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
#include <google/protobuf/message.h>
#include <google/protobuf/text_format.h>
#include "utils/exception.hpp"
#include "proto/utils.hpp"

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

google::protobuf::FieldDescriptor const*
get_oneof_field_desc(const google::protobuf::Message& msg,
                     const std::string& oneof_name)
{
  auto desc = msg.GetDescriptor();
  auto oneof_handle = desc->FindOneofByName(oneof_name);
  if (oneof_handle == nullptr)
  {
    std::string err_str
      = "Unable to identify the type 'oneof " + oneof_name + "' in message {"
      + desc->DebugString() + "}\n";
    WCS_THROW(err_str);
  }
  auto reflex = msg.GetReflection();

  return reflex->GetOneofFieldDescriptor(msg, oneof_handle);
}

bool has_oneof(google::protobuf::Message const& msg,
               std::string const& oneof_name)
{
  return (get_oneof_field_desc(msg, oneof_name) != nullptr);
}

const google::protobuf::Message&
get_oneof_message(const google::protobuf::Message& msg,
                  const std::string& oneof_name)
{
  auto oneof_field = get_oneof_field_desc(msg, oneof_name);
  if (oneof_field == nullptr) {
    std::string err_str = "The value of 'oneof " + oneof_name
                        + "' has not been set in message {"
                        + msg.DebugString() + "}\n";
    WCS_THROW(err_str);
  }

  if (oneof_field->type() != google::protobuf::FieldDescriptor::TYPE_MESSAGE) {
    WCS_THROW("Oneof field is not of message type.");
  }

  return msg.GetReflection()->GetMessage(msg, oneof_field);
}

} // end of namespace wcs
