// Whole Cell Model Simulator
/** @file
 * \ingroup wcs_reaction_network
 * \brief vertex definition for the reaction graph representation
 */
/** \ingroup wcs_reaction_network
 * \class wcs::VertexPropertyBase
 * \brief vertex definition for the reaction graph representation
 *
 * The bundled vertex property data structure used to define the BGL graph,
 * which represents a reaction network.
 */
#ifndef __WCS_VERTEX_PROPERTY_BASE_HPP__
#define __WCS_VERTEX_PROPERTY_BASE_HPP__

#include "wcs_types.hpp"
#include <iostream>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

class VertexPropertyBase {
 public:
  VertexPropertyBase() {} ///< BGL requires it to be default constructable
  VertexPropertyBase(const VertexPropertyBase& rhs) = default;
  VertexPropertyBase(VertexPropertyBase&& rhs) noexcept;
  VertexPropertyBase& operator=(const VertexPropertyBase& rhs) = default;
  VertexPropertyBase& operator=(VertexPropertyBase&& rhs) noexcept;
  virtual ~VertexPropertyBase() = 0;
  std::unique_ptr<VertexPropertyBase> clone() const;

 protected:
  void reset(VertexPropertyBase& obj);

 private:
  virtual VertexPropertyBase* clone_impl() const = 0;

 friend ::wcs::GraphFactory;

 template <typename G>
 friend std::ostream& ::wcs::write_graphviz_of_any_vertex_list(std::ostream& os, const G& g);
};

/**@}*/
} // end of namespace wcs
#endif // __WCS_VERTEX_PROPERTY_BASE_HPP__
