include config.mk

.PHONY: all target_dep target_rule test clean

#TARGET := reaction
TARGET := ssa
SRCS := $(TARGET).cpp

WCM_SUBDIRS = reaction_network utils sim_methods
WCM_OBJS = reaction_network/*.o utils/graph_factory.o utils/exception.o utils/trace.o utils/trace_ssa.o sim_methods/ssa_nrm.o

all: target_dep target_rule
	@for i in ${WCM_SUBDIRS}; \
	do \
	    if [ -f $$i/Makefile ]; \
	    then \
	        echo "Making $$i ..."; \
	        ($(MAKE) -f $$i/Makefile all) || exit 1; \
	    fi; \
	done
	@rm -f $(TARGET)
	$(MAKE) $(TARGET)


        #(cd $$i && $(MAKE) all) || exit 1; \


$(TARGET): $(SRCS) $(WCM_OBJS)
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -o $@ $@.cpp $(WCM_OBJS) $(LIBS)


%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

%.d: %.cpp
	@echo "$@ `$(CXXDEP) $(CXXFLAGS) $(CPPFLAGS) -MM $<`" > $@

target_dep: $(addsuffix .d, $(basename $(SRCS)))

target_rule: $(OBJS)

test: ${TARGET}
	./$(TARGET) utils/unit_tests/examples/simple.graphml -g simple_poly.dot -o result.txt -i 20 -s 47
	dot -Teps simple_poly.dot -o simple_poly.eps
#	./${TARGET} > graph.dot
#	dot -Teps graph.dot -o graph.eps

clean:
	@rm -f $(TARGET)
	@for i in ${WCM_SUBDIRS}; \
	 do \
	   if [ -f $$i/Makefile ]; \
	   then \
	     echo "Cleaning $$i ..."; \
	     ($(MAKE) -f $$i/Makefile clean) || exit 1; \
	   fi; \
	 done
	@rm -f *.core
	@rm -f *.d *.o
	@rm -f *.dot *.eps
