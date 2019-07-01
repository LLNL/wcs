

set(CLOCK gtod CACHE STRING "Clock setting to use, gtod or amd64")
set(ROSS_DEFINES "-D_GNU_SOURCE" CACHE STRING "Extra defines for ROSS compilation")

blt_add_library(
  NAME ROSS
  SOURCES
    queue/splay.cc
    queue/splay.cc
    ross-random.cc
    rand-clcg4.cc
    clock/${CLOCK}.cc
    gvt/mpi_allreduce.cc
    network-mpi.cc
    hash-quadratic.cc
    buddy.cc
    lz4.cc
    tw-state.cc
    ross-global.cc
    tw-event.cc
    tw-opts.cc
    tw-pe.cc
    tw-kp.cc
    tw-lp.cc
    tw-timing.cc
    tw-sched.cc
    tw-setup.cc
    tw-signal.cc
    tw-stats.cc
    tw-util.cc
    avl_tree.cc
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/../include
  DEFINES ${ROSS_DEFINES}
  DEPENDS_ON mpi
)
