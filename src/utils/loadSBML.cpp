#if defined(WCS_HAS_CONFIG)
#include "wcs_config.hpp"
#else
#error "no config"
#endif

#if defined(WCS_HAS_SBML)
#include <iostream>

#include <sbml/SBMLTypes.h>
#include <sbml/common/extern.h>


LIBSBML_CPP_NAMESPACE_USE

BEGIN_C_DECLS

int main(int argc, char** argv)
{
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " filename" << std::endl;
    return 0;
  }

  const char* filename = argv[1];
  SBMLReader reader;
  SBMLDocument* document = reader.readSBML(filename);

  const unsigned num_errors = document->getNumErrors();

  if (num_errors > 0u) {
    std::cout << num_errors << " error(s) in reading "
              << filename << std::endl;

    document->printErrors(std::cerr);

    delete document;
    return static_cast<int>(num_errors);
  } else {
    std::cout << "SBML file " << filename << ":" << std::endl
              << "\tLevel = " << document->getLevel()
              << "\tVersion = " << document->getVersion()
              << std::endl;
  }

  Model* model = document->getModel();

  if (model == NULL) {
    std::cout << "Faile to get model from " << filename << std::endl;
    delete document;
    return -1;
  }

  if (model->isSetSBOTerm()) {
    std::cout << "Num species: " << model->getNumSpecies() << std::endl;
    std::cout << "Num reactions: " << model->getNumReactions() << std::endl;
  }

  delete document;
  return 0;
}

END_C_DECLS
#endif // defined(WCS_HAS_SBML)
