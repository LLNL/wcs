#include "reaction_network/species.hpp"
#include <limits>

namespace wcs {
/** \addtogroup wcs_reaction_network
 *  *  @{ */

species_cnt_t Species::m_max_count = std::numeric_limits<species_cnt_t>::max();

Species::Species()
: VertexPropertyBase(),
  m_count(static_cast<species_cnt_t>(0))
{}

Species::Species(const Species& rhs)
: VertexPropertyBase(rhs),
  m_count(rhs.m_count)
{}

Species::Species(Species&& rhs) noexcept
: VertexPropertyBase(std::move(rhs)),
  m_count(rhs.m_count)
{
  if (this != &rhs) {
    reset(rhs);
  }
}

Species& Species::operator=(const Species& rhs)
{
  if (this != &rhs) {
    VertexPropertyBase::operator=(rhs);
    m_count = rhs.m_count;
  }
  return *this;
}

Species& Species::operator=(Species&& rhs) noexcept
{
  if (this != &rhs) {
    VertexPropertyBase::operator=(std::move(rhs));
    m_count = rhs.m_count;
    reset(rhs);
  }
  return *this;
}

Species::~Species() {}

std::unique_ptr<Species> Species::clone() const
{
  return std::unique_ptr<Species>(this->clone_impl());
}

Species* Species::clone_impl() const
{
  return (new Species(*this));
}

void Species::reset(Species& obj)
{
  VertexPropertyBase::reset(obj);
  obj.m_count = static_cast<species_cnt_t>(0);
}

bool Species::inc_count()
{
  if (m_count >= m_max_count) {
    return false;
  }
  m_count ++;
  return true;
}

bool Species::dec_count()
{
  if (m_count <= 0) {
    return false;
  }
  m_count --;
  return true;
}

bool Species::inc_count(const species_cnt_t c)
{
  if ((m_max_count -  m_count) < c) {
    return false;
  }
  m_count += c;
  return true;
}

bool Species::dec_count(const species_cnt_t c)
{
  if (m_count < c) {
    return false;
  }
  m_count -= c;
  return true;
}

bool Species::set_count(const species_cnt_t c)
{
  if ((m_count < 0) || (m_count > m_max_count)) {
    return false;
  }
  m_count = c;
  return true;
}

species_cnt_t Species::get_count() const
{
  return m_count;
}

/**@}*/
} // end of namespace wcs
