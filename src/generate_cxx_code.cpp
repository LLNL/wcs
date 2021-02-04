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

#include <unordered_map>
#include <cstdio>

#include "utils/generate_cxx_code.hpp"

#if defined(WCS_HAS_SBML)
#include <iostream>

#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>


LIBSBML_CPP_NAMESPACE_USE


int main(int argc, char** argv)
{
  if ((argc != 2) && (argc != 3))  {
    std::cout << "Usage: " << argv[0] << " filename [gen_library(0|1)]" << std::endl;
    return 0;
  }

  const char* filename = argv[1];
  const bool gen_lib = (argc > 2) && (atoi(argv[2]) != 0);
  SBMLReader reader;
  SBMLDocument* document = reader.readSBML(filename);
  const unsigned int num_errors = document->getNumErrors();

  if (num_errors > 0u) {
    std::cout << num_errors << " error(s) in reading "
              << filename << std::endl;

    document->printErrors(std::cerr);

    delete document;
    return static_cast<int>(num_errors);
  }

  const Model* model = document->getModel();

  if (model == nullptr) {
    std::cout << "Failed to get model from " << filename << std::endl;
    delete document;
    return -1;
  }
  using params_map = std::unordered_map <std::string, std::vector<std::string>>;
  using rate_rules_dep = std::unordered_map <std::string, std::set<std::string>>;
  params_map dep_params_f, dep_params_nf;
  rate_rules_dep rate_rules_dep_map;
  wcs::generate_cxx_code code_generator("wcs_generated_code.so", true, false, false);
  code_generator.generate_code(*model, dep_params_f, dep_params_nf,
                               rate_rules_dep_map);
  std::cout << "Generated filename: " << code_generator.get_src_filename() << "\n";
  if (gen_lib) {
    const std::string library_file = code_generator.compile_code();
    std::cout << "library file: " << library_file << std::endl;
  }

  delete document;
  return 0;
}
#endif // defined(WCS_HAS_SBML)
