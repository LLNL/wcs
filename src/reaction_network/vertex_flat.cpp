/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "reaction_network/vertex_flat.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

std::map<VertexFlat::vertex_type, std::string> VertexFlat::vt_str {
  {VertexFlat::_undefined_, "undefined"},
  {VertexFlat::_species_,   "species"},
  {VertexFlat::_reaction_,  "reaction"}
};

std::map<std::string, VertexFlat::vertex_type> VertexFlat::str_vt {
  {"undefined", VertexFlat::_undefined_},
  {"species",   VertexFlat::_species_},
  {"reaction",  VertexFlat::_reaction_},
  {"",          VertexFlat::_undefined_}
};

VertexFlat::VertexFlat()
: m_type(_undefined_),
  m_label(""),
  m_pid(unassigned_partition),
  m_count(static_cast<s_cnt_t>(0)),
  m_rate(static_cast<r_rate_t>(0)),
  m_rate_formula("")
{
  m_typeid = static_cast<int>(m_type);
}

/*
VertexFlat::VertexFlat(const vertex_type vt, const std::string& lb)
: m_type(vt), m_label(lb),
  m_pid(unassigned_partition),
  m_count(static_cast<s_cnt_t>(0)),
  m_rate(static_cast<r_rate_t>(0)),
  m_rate_formula("")
{
  m_typeid = static_cast<int>(m_type);
}
*/

void VertexFlat::set_type(const vertex_type vt)
{
  m_type = vt;
  m_typeid = static_cast<int>(m_type);
}

/**
 * Convert the vertex type id in integer format, which is read from GraphML
 * input, into the enum type that is internally used.
 */
void VertexFlat::set_type()
{
  m_type = static_cast<vertex_type>(m_typeid);
}

VertexFlat::vertex_type VertexFlat::get_type() const
{
  return m_type;
}

int VertexFlat::get_typeid() const
{
  return m_typeid;
}

std::string VertexFlat::get_type_str() const
{
  return vt_str[m_type];
}

void VertexFlat::set_label(const std::string& lb)
{
  m_label = lb;
}

std::string VertexFlat::get_label() const
{
  return m_label;
}

void VertexFlat::set_partition(const partition_id_t pid)
{
  m_pid = pid;
}

partition_id_t VertexFlat::get_partition() const
{
  return m_pid;
}

bool VertexFlat::inc_count()
{
  m_count ++;
  return true;
}

bool VertexFlat::dec_count()
{
  if (m_count <= 0) {
    return false;
  }
  m_count --;
  return true;
}

bool VertexFlat::set_count(VertexFlat::s_cnt_t c)
{
  if (m_count < 0) {
    return false;
  }
  m_count = c;
  return true;
}

VertexFlat::s_cnt_t VertexFlat::get_count() const
{
  return m_count;
}

void VertexFlat::set_rate_constant(VertexFlat::r_rate_t k)
{
  m_rate_const = k;
}

VertexFlat::r_rate_t VertexFlat::get_rate_constant() const
{
  return m_rate_const;
}

VertexFlat::r_rate_t VertexFlat::compute_rate(const std::vector<VertexFlat::s_cnt_t>& cnts)
{
  s_cnt_t pop_factor = static_cast<s_cnt_t>(1);
  for(const auto c : cnts) {
    pop_factor *= c;
  }
  m_rate = m_rate_const * pop_factor;

  return m_rate;
}

VertexFlat::r_rate_t VertexFlat::get_rate() const
{
  return m_rate;
}

void VertexFlat::set_rate_formula(const std::string& f)
{
  m_rate_formula = f;
}

const std::string& VertexFlat::get_rate_formula() const
{
  return m_rate_formula;
}

/**@}*/
} // end of namespace wcs
