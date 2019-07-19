#OPT += -Wno-unknown-pragmas
include config.mk


SRCS := reaction_network/vertex_flat.cpp \
        reaction_network/vertex.cpp \
        reaction_network/edge.cpp \
        reaction_network/vertex_property_base.cpp \
        reaction_network/reaction_base.cpp \
        reaction_network/species.cpp \
        reaction_network/network.cpp

OBJS := $(addsuffix .o, $(basename $(SRCS)))

.PHONY: all clean target_dep target_rule

all: target_dep target_rule


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

%.d: %.cpp
	@echo "$@ `$(CXXDEP) $(CXXFLAGS) $(CPPFLAGS) -MM $<`" > $@

clean:
	@rm -f reaction_network/*.o reaction_network/*.d

target_dep: $(addsuffix .d, $(basename $(SRCS)))

target_rule: $(OBJS)
	
ifneq ($(MAKECMDGOALS),clean)
    -include $(OBJS:.o=.d)
endif
