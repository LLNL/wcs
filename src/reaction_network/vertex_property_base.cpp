/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "reaction_network/vertex_property_base.hpp"
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

VertexPropertyBase::VertexPropertyBase(VertexPropertyBase&& rhs) noexcept
{
  if (this != &rhs) {
  }
}

VertexPropertyBase& VertexPropertyBase::operator=(VertexPropertyBase&& rhs) noexcept
{
  if (this != &rhs) {
    reset(rhs);
  }
  return *this;
}

VertexPropertyBase::~VertexPropertyBase()
{}

std::unique_ptr<VertexPropertyBase> VertexPropertyBase::clone() const
{
  return std::unique_ptr<VertexPropertyBase>(this->clone_impl());
}

VertexPropertyBase* VertexPropertyBase::clone_impl() const
{
  WCS_THROW("Method not implemented in abstract class.");
  return nullptr;//(new VertexPropertyBase(*this));
}

void VertexPropertyBase::reset(VertexPropertyBase& obj)
{}

/**@}*/
} // end of namespace wcs
