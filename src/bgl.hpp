#ifndef __WCS_BGL_HPP__
#define __WCS_BGL_HPP__

#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
// To suppress the gcc compiler warning 'maybe-uninitialized'
// from the boost graph source code.
// clang does not recognize this particular diagnostic flag.
#include <boost/graph/adjacency_list.hpp>
#if !defined(__clang__) && defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#endif // __WCS_BGL_HPP__
