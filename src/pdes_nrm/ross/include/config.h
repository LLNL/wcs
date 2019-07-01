#ifndef ROSS_CONFIG__
#define ROSS_CONFIG__

// ROSS Configuration Options
#define ROSS_QUEUE_splay
#define ROSS_RAND_clcg4
#define ROSS_NETWORK_mpi
#define ROSS_GVT_mpi_allreduce
#define ROSS_CLOCK_gtod
#define ARCH_i386

// ROSS Core
#define HAVE_CTIME 1
//#define ROSS_VERSION "${GIT_SHA1}"
#define AVL_TREE 1
#define AVL_NODE_COUNT (1<<18)

//#cmakedefine USE_BGPM
//#cmakedefine ROSS_MEMORY
#define RAND_NORMAL
#define ROSS_timing
#define ROSS_runtime_checks

#endif
