/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef	 __WCS_SIM_METHODS_UPDATE_HPP__
#define	 __WCS_SIM_METHODS_UPDATE_HPP__
#include <vector>
#include "wcs_types.hpp"
#include "reaction_network/network.hpp"


namespace wcs {
/** \addtogroup wcs_utils
 *  @{ */

/// Type of the BGL vertex descriptor for the species vertex
using s_desc_t = wcs::Network::v_desc_t;

/** Type for the update based on the species count
 *  To reduce memory footprint, the update can be of a short signed intgeter
 *  type.
 */
using cnt_update_t = std::pair<s_desc_t, wcs::stoic_t>;
// If stoic_t is an unsigned type, use species_cnt_diff_t instead
//using cnt_update_t = std::pair<s_desc_t, wcs::species_cnt_diff_t>;

/// Type for the update based on the species concentration
using conc_update_t = std::pair<s_desc_t, wcs::concentration_t>;


using cnt_updates_t = std::vector<cnt_update_t>;
using conc_updates_t = std::vector<conc_update_t>;

/**@}*/
} // end of namespace wcs
#endif // __WCS_SIM_METHODS_UPDATE_HPP__
