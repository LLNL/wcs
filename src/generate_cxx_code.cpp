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

#include "utils/file.hpp"
#include "utils/generate_cxx_code.hpp"

#if defined(WCS_HAS_SBML)
#include <iostream>

#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>


LIBSBML_CPP_NAMESPACE_USE


int main(int argc, char** argv)
{
  if ((argc < 2) || (argc > 6))  {
    std::cout << "Usage: " << argv[0]
              << " model_filename [gen_library(0|1)"
              << " [compilation_error_log(0|1) [chunk_size [tmpdir]]]]" << std::endl;
    return EXIT_SUCCESS;
  }

  const char* model_filename = argv[1];
  const bool gen_lib = (argc > 2) && (atoi(argv[2]) != 0);
  const bool show_error = (argc > 3) && (atoi(argv[3]) != 0);
  const unsigned int chunk_size = ((argc > 4)? atoi(argv[4]) : 1000u);
  const std::string tmp_dir = ((argc > 5)? argv[5] : "/tmp");
  SBMLReader reader;
  SBMLDocument* document = reader.readSBML(model_filename);
  const unsigned int num_errors = document->getNumErrors();

  if (num_errors > 0u) {
    std::cout << num_errors << " error(s) in reading "
              << model_filename << std::endl;

    document->printErrors(std::cerr);

    delete document;
    return static_cast<int>(num_errors);
  }

  const Model* model = document->getModel();

  if (model == nullptr) {
    std::cout << "Failed to get model from " << model_filename << std::endl;
    delete document;
    return EXIT_FAILURE;
  }

  std::cout << "chunk size: " << chunk_size << std::endl;

  const std::string lib_filename = wcs::get_libname_from_model(model_filename);
  wcs::generate_cxx_code code_generator(lib_filename, true, show_error, false,
                                        tmp_dir, chunk_size);

  using params_map = std::unordered_map <std::string, std::vector<std::string>>;
  using rate_rules_dep = std::unordered_map <std::string, std::set<std::string>>;
  params_map dep_params_f, dep_params_nf;
  rate_rules_dep rate_rules_dep_map;

  code_generator.generate_code(*model, dep_params_f, dep_params_nf,
                               rate_rules_dep_map);

  std::cout << "Generated source filenames:";
  for (const auto& fn: code_generator.get_src_filenames()) {
    std::cout << ' ' << fn;
  }
  std::cout  << std::endl;

  if (gen_lib) {
    code_generator.compile_code();
    std::cout << "library file: " << lib_filename << std::endl;
  }

  delete document;
  return EXIT_SUCCESS;
}
#endif // defined(WCS_HAS_SBML)
