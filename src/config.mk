# works with clang/6.0.0
#BOOST_ROOT = /usr/tce/packages/boost/boost-1.69.0-mvapich2-2.2-gcc-4.9.3
# works with gcc/7.1.0
#BOOST_ROOT = /usr/tce/packages/boost/boost-1.66.0-mvapich2-2.2-gcc-6.1.0

# The cmake recipe is under development. Until that becomes available,
# make sure to export CXX, BOOST_ROOT and 
# LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/${BOOST_ROOT}/lib
# Then, type make

CPPFLAGS += -I. -I${BOOST_ROOT}/include
CXXFLAGS += -Wall -O3 -std=c++17
LIBS += -lstdc++ -L${BOOST_ROOT}/lib -lboost_graph -lboost_regex
LIBS += -lboost_filesystem -lboost_system
CXXDEP := $(CXX) -E
