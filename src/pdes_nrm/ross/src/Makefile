CC= mpicxx

# Compiler flags
CFLAGS= -I../include
OPT= -g -O2 -Wall


SHELL= /bin/sh

# Target library name
TARGET= libROSS.a



# IBM/BGQ:
ARCH= i386
CLOCK= amd64
QUEUE= splay

# arm/Linux-on-Android:
ARCH= i386
CLOCK= gtod
QUEUE= splay

#ROSS_FLAGS= -D_GNU_SOURCE -DROSS_CLOCK_$(CLOCK) -DARCH_$(ARCH)
ROSS_FLAGS= -D_GNU_SOURCE
# -DROSS_OPTION_LIST='{}'

SRCS= \
  queue/splay.cc  ross-random.cc        rand-clcg4.cc   clock/${CLOCK}.cc \
  gvt/mpi_allreduce.cc      network-mpi.cc hash-quadratic.cc buddy.cc          \
  lz4.cc             tw-state.cc           ross-global.cc    tw-event.cc       \
  tw-opts.cc         tw-pe.cc              tw-kp.cc          tw-lp.cc          \
  tw-timing.cc       tw-sched.cc           tw-setup.cc       tw-signal.cc      \
  tw-stats.cc        tw-util.cc            avl_tree.cc


# Create names of object files
OBJS= $(addsuffix .o, $(basename $(SRCS)))

# Rule to compile C source file into object file
%.o: %.cc
	$(CC) $(CFLAGS) $(ROSS_FLAGS) $(OPT) -o $@ -c $<
# Create dependency file from C source
%.d: %.cc
	(echo -n "$@ " ; $(CC) $(CFLAGS) $(ROSS_FLAGS) -MM $< ) > $@

# Default rule
all: target_dep target_rule

# Clean objects, dependency files, and target library
clean:
	rm -f $(OBJS) $(OBJS:.o=.d) $(TARGET)

# Install (copy) target library into ../lib/
install: $(TARGET)
	mkdir -p ../lib
	cp -a $(TARGET) ../lib/

# Create names of depenency files
target_dep: $(addsuffix .d, $(basename $(SRCS)))

# Link library
target_rule: $(OBJS)
	rm -rf $(TARGET)
	ar rc $(TARGET) $(OBJS)

# Include dependency files
ifneq ($(MAKECMDGOALS),clean)
-include $(OBJS:.o=.d)
endif
