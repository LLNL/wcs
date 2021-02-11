/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef __USE_GNU // for CPU_ISSET
#define __USE_GNU 1
#endif

#include <cstdio>   // perror
#include <unistd.h> // sysconf
#include <sched.h>  // sched_getaffinity
#include <iostream>
#include <unordered_map>
#include <omp.h>
#include "utils/omp_diagnostics.hpp"


namespace wcs {

std::string get_omp_version()
{
#if defined(_OPENMP)
  static const std::unordered_map<unsigned, std::string> omp_versions{
    {200505,"2.5"},{200805,"3.0"},{201107,"3.1"},{201307,"4.0"},
    {201511,"4.5"},{201811,"5.0"},{202011,"5.1"},{201811,"5.0"}};
  const auto it = omp_versions.find(_OPENMP);
  if (it == omp_versions.cend()) {
    return std::to_string(_OPENMP);
  }
  return std::to_string(_OPENMP) + " " + it->second;
#else
  return "NA";
#endif // defined(_OPENMP)
}

// Adopted from OpenMP Diagnostic code from Edgar Leon at LLNL
/// Get number of processing units (cores or hardware threads)
static int get_num_pus()
{
  int pus = 0;
  if ((pus = sysconf(_SC_NPROCESSORS_ONLN)) < 0) {
    perror("sysconf");
  }
  return pus;
}

/// Get the affinity
static int get_affinity(std::vector<uint8_t>& cpus)
{
  int rc = 0;
  cpu_set_t resmask;
  int pus = get_num_pus();

  CPU_ZERO(&resmask);
  if ((rc = sched_getaffinity(0, sizeof(resmask), &resmask)) < 0) {
    perror("sched_getaffinity");
    return rc;
  }

  for (int i = 0; i < pus; i++) {
    if (CPU_ISSET(i, &resmask)) {
      cpus.push_back(static_cast<uint8_t>(i));
    }
  }

  return rc;
}

void my_omp_affinity::get()
{
  m_ancestor_id.clear();
#if defined(_OPENMP)
  m_tid = omp_get_thread_num();
  m_num_threads = omp_get_num_threads();
  m_my_level = omp_get_level();
  m_ancestor_id.resize(m_my_level-1);
  for (int i = 1; i < m_my_level; i++) {
    m_ancestor_id[i] = omp_get_ancestor_thread_num(i);
  }
#else
  m_tid = 1;
  m_num_threads = 1;
#endif

  m_cpus.reserve(get_num_pus());
  get_affinity(m_cpus);
}

void my_omp_affinity::print() const
{
  std::string id_str = "OMP thread ";

  if (m_my_level > 1) {
    id_str += "at level " + std::to_string(m_my_level) + " ";
    for (int i = 1; i < m_my_level; i++) {
      id_str += std::string((m_ancestor_id[i] < 10)? "0" : "")
                 + std::to_string(m_ancestor_id[i]) + "-";
    }
  }
  id_str += std::string((m_tid < 10)? "0" : "") + std::to_string(m_tid);

  std::string msg
    = id_str + "/" + std::to_string(m_num_threads) + " can run on cpu ";

  for (size_t i = 0u; i + 1 < m_cpus.size(); i++) {
    msg += std::to_string(m_cpus[i]) + ", ";
  }
  if (m_cpus.size() > 0u) {
    msg += std::to_string(m_cpus.back()) + "\n";
  }
#if defined(_OPENMP)
  if (omp_get_level()) {
    #pragma omp critical
    {
      std::cout << msg;
    }
  } else {
    std::cout << msg;
  }
#else
  std::cout << msg;
#endif // defined(_OPENMP)
}

std::string to_string_omp_schedule_kind(int kind)
{
#if defined(_OPENMP)
  switch (static_cast<unsigned>(kind))
  {
    case static_cast<unsigned>(0x1): return "omp_sched_static";
    case static_cast<unsigned>(0x2): return "omp_sched_dynamic";
    case static_cast<unsigned>(0x3): return "omp_sched_guided";
    case static_cast<unsigned>(0x4): return "omp_sched_auto";
    case static_cast<unsigned>(0x80000000u): return "omp_sched_auto";
  }
#endif // defined(_OPENMP)
  return "unknown";
}

void set_static_schedule()
{
#if defined(_OPENMP)
  omp_sched_t kind = omp_sched_static;
  int chunk_size = 1;
  omp_get_schedule(&kind, &chunk_size);

  std::string msg;
  if (kind != omp_sched_static) {
    omp_set_schedule(omp_sched_static, 1);
    msg = "change " + to_string_omp_schedule_kind(kind)
        + " to omp_sched_static";
  } else {
    msg = "omp_sched_static with chunk size " + std::to_string(chunk_size);
  }
  std::cout << msg << std::endl;
#endif // defined(_OPENMP)
}

} // namespace wcs
