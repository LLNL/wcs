/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "reaction_network/vertex.hpp"
#include "reaction_network/edge.hpp"
#include <memory>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

std::map<Vertex::vertex_type, std::string> Vertex::vt_str {
  {Vertex::_undefined_, "undefined"},
  {Vertex::_species_,   "species"},
  {Vertex::_reaction_,  "reaction"}
};

std::map<std::string, Vertex::vertex_type> Vertex::str_vt {
  {"undefined", Vertex::_undefined_},
  {"species",   Vertex::_species_},
  {"reaction",  Vertex::_reaction_},
  {"",          Vertex::_undefined_}
};

Vertex::Vertex()
: m_type(_undefined_), m_label(""),
  m_p(nullptr)
{
  m_typeid = static_cast<int>(m_type);
}

Vertex::Vertex(const Vertex& rhs)
: m_type(rhs.m_type),
  m_typeid(rhs.m_typeid),
  m_label(rhs.m_label),
  m_p((!rhs.m_p)? nullptr : rhs.m_p.get()->clone())
{
}

Vertex::Vertex(Vertex&& rhs) noexcept
: m_type(rhs.m_type),
  m_typeid(rhs.m_typeid)
{
  if (this != &rhs) {
    m_label = std::move(rhs.m_label);
    m_p = std::move(rhs.m_p);

    reset(rhs);
  }
}

Vertex& Vertex::operator=(const Vertex& rhs)
{
  if (this != &rhs) {
    m_type = rhs.m_type;
    m_typeid = rhs.m_typeid;
    m_label = rhs.m_label;
    m_p = (!rhs.m_p)? nullptr : rhs.m_p.get()->clone();
  }
  return *this;
}

Vertex& Vertex::operator=(Vertex&& rhs) noexcept
{
  if (this != &rhs) {
    m_type = rhs.m_type;
    m_typeid = rhs.m_typeid;
    m_label = std::move(rhs.m_label);
    m_p = std::move(rhs.m_p);

    reset(rhs);
  }
  return *this;
}

Vertex::~Vertex()
{}

std::unique_ptr<Vertex> Vertex::clone() const
{
  return std::unique_ptr<Vertex> (this->clone_impl());
}

Vertex* Vertex::clone_impl() const
{
  return (new Vertex(*this));
}

void Vertex::reset(Vertex& obj)
{
   obj.m_label.clear();
   obj.m_p = nullptr;
}

void Vertex::set_type(const Vertex::vertex_type vt)
{
  m_type = vt;
  m_typeid = static_cast<int>(m_type);
}

void Vertex::set_type()
{
  m_type = static_cast<vertex_type>(m_typeid);
}

Vertex::vertex_type Vertex::get_type() const
{
  return m_type;
}

int Vertex::get_typeid() const
{
  return m_typeid;
}

std::string Vertex::get_type_str() const
{
  return vt_str[m_type];
}

void Vertex::set_label(const std::string& lb)
{
  m_label = lb;
}

std::string Vertex::get_label() const
{
  return m_label;
}


/**@}*/
} // end of namespace wcs
