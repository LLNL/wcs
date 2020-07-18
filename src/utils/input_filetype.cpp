/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "utils/input_filetype.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>


namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

input_filetype::input_filetype(const std::string filename)
: m_filename(filename)
{}

input_filetype::input_type input_filetype::detect() const
{
  std::ifstream file;
  file.open(m_filename);
  std::string line;
  const std::string commentline("<!--");
  const std::string graphmlline("graphml");
  const std::string sbmlline("sbml");

  if (!file) //checks to see if file opens properly
  {
    file.close(); // Remember to close the file.
    return input_filetype::input_type::_ioerror_;
  }
  else
  {
    for(int i=0; i<10; i++) {
      if (std::getline(file, line)){
        std::transform(line.begin(), line.end(), line.begin(),
        [](unsigned char c){ return std::tolower(c); });
        const size_t pos = line.find(commentline);
        const size_t pos1 = line.find(graphmlline);
        const size_t pos2 = line.find(sbmlline);
        if (pos != std::string::npos) {
          // If there is a comment line, it decreases the i in order to
          // check only the first 10 lines of a SBML or a graphML file
          i--;
        } else if (pos1 != std::string::npos) {  ///graphml file
          file.close(); // Remember to close the file.
          return input_filetype::input_type::_graphml_;
        } else if (pos2 != std::string::npos) {  ///sbml file
          file.close(); // Remember to close the file.
          return input_filetype::input_type::_sbml_;
        }
      }
    }

  }
  file.close(); // Remember to close the file.
  return input_filetype::input_type::_unknown_;
}


/**@}*/
} // end of namespace wcs
