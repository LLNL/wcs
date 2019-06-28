#OPT += -Wno-unknown-pragmas
include config.mk


SRCS := sim_methods/ssa_nrm.o

OBJS := $(addsuffix .o, $(basename $(SRCS)))

.PHONY: all clean target_dep target_rule

all: target_dep target_rule


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

%.d: %.cpp
	@echo "$@ `$(CXXDEP) $(CXXFLAGS) $(CPPFLAGS) -MM $<`" > $@

clean: 
	@rm -f sim_methods/*.o sim_methods/*.d

target_dep: $(addsuffix .d, $(basename $(SRCS)))

target_rule: $(OBJS)
	
ifneq ($(MAKECMDGOALS),clean)
    -include $(OBJS:.o=.d)
endif
