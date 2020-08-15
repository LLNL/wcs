/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#include "reaction_network/reaction_base.hpp"
#include "utils/exception.hpp"

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  @{ */

ReactionBase::ReactionBase()
: VertexPropertyBase(),
  m_rate(static_cast<reaction_rate_t>(0)),
  m_rate_const(static_cast<reaction_rate_t>(0)),
  m_rate_formula("")
{
  set_calc_rate_fn();
}

ReactionBase::ReactionBase(const ReactionBase& rhs)
: VertexPropertyBase(rhs),
  m_rate(rhs.m_rate),
  m_rate_const(rhs.m_rate_const),
  m_rate_formula(rhs.m_rate_formula),
  m_calc_rate(rhs.m_calc_rate)
{}

ReactionBase::ReactionBase(ReactionBase&& rhs) noexcept
: VertexPropertyBase(std::move(rhs)),
  m_rate(rhs.m_rate),
  m_rate_const(rhs.m_rate_const)
{
  if (this != &rhs) {
    m_rate_formula = std::move(rhs.m_rate_formula);
    m_calc_rate = std::move(rhs.m_calc_rate);

    reset(rhs);
  }
}

ReactionBase& ReactionBase::operator=(const ReactionBase& rhs)
{
  if (this != &rhs) {
    VertexPropertyBase::operator=(rhs);
    m_rate = rhs.m_rate;
    m_rate_const = rhs.m_rate_const;
    m_rate_formula = rhs.m_rate_formula;
    m_calc_rate = rhs.m_calc_rate;
  }
  return *this;
}

ReactionBase& ReactionBase::operator=(ReactionBase&& rhs) noexcept
{
  if (this != &rhs) {
    VertexPropertyBase::operator=(std::move(rhs));
    m_rate = rhs.m_rate;
    m_rate_const = rhs.m_rate_const;
    m_rate_formula = std::move(rhs.m_rate_formula);
    m_calc_rate = std::move(rhs.m_calc_rate);

    reset(rhs);
  }
  return *this;
}

ReactionBase::~ReactionBase(){}

std::unique_ptr<ReactionBase> ReactionBase::clone() const
{
  WCS_THROW("Method not implemented in abstract class.");
  return std::unique_ptr<ReactionBase>(nullptr);
  //return std::unique_ptr<ReactionBase>(this->clone_impl());
}

/*
ReactionBase* ReactionBase::clone_impl() const
{
  return (new ReactionBase(*this));
}
*/

void ReactionBase::reset(ReactionBase& obj)
{
  obj.m_rate = static_cast<reaction_rate_t>(0);
  obj.m_rate_const = static_cast<reaction_rate_t>(0);
  obj.m_rate_formula.clear();
  obj.m_calc_rate = nullptr;
}

void ReactionBase::ReactionBase::set_rate_constant(reaction_rate_t k)
{
  m_rate_const = k;
}

reaction_rate_t ReactionBase::get_rate_constant() const
{
  return m_rate_const;
}

void ReactionBase::set_calc_rate_fn()
{
  m_calc_rate = [] (const std::vector<reaction_rate_t>& params)
    {
      reaction_rate_t rate = static_cast<reaction_rate_t>(1);
      for(const auto p : params) {
        rate *= p;
      }
      return rate;
    };
}

void ReactionBase::set_calc_rate_fn(
  const std::function<
    reaction_rate_t (const std::vector<reaction_rate_t>&)
  >& calc_rate)
{
  m_calc_rate = calc_rate;
}

reaction_rate_t ReactionBase::calc_rate(std::vector<reaction_rate_t>&& params)
{
  params.push_back(m_rate_const);
  m_rate = (!m_calc_rate)? 0.0 : m_calc_rate(params);
  return m_rate;
}

void ReactionBase::set_rate(const reaction_rate_t rate)
{
  m_rate = rate;
}

reaction_rate_t ReactionBase::get_rate() const
{
  return m_rate;
}

void ReactionBase::set_rate_formula(const std::string& f)
{
  m_rate_formula = f;
}

const std::string& ReactionBase::get_rate_formula() const
{
  return m_rate_formula;
}

/**@}*/
} // end of namespace wcs
