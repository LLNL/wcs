BASE_DIR = ../..
CUR_DIR := $(BASE_DIR)/utils/unit_tests
include $(BASE_DIR)/config.mk
CPPFLAGS += -I$(BASE_DIR)
LIBS += -lboost_filesystem -lboost_system

.PHONY: all target_dep target_rule test clean

TARGET = t_graph_factory
#TARGET = t_rngen
SRCS = $(TARGET).cpp

WCM_OBJS = $(BASE_DIR)/reaction_network/*.o $(BASE_DIR)/utils/graph_factory.o $(BASE_DIR)/utils/exception.o

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

#	(cd $$i && $(MAKE) all) || exit 1; \

${TARGET}: $(SRCS) $(WCM_OBJS)
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -o $@ $@.cpp $(WCM_OBJS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

%.d: %.cpp
	@echo "$@ `$(CXXDEP) $(CXXFLAGS) $(CPPFLAGS) -MM $<`" > $@

target_dep: $(addsuffix .d, $(basename $(SRCS)))

target_rule: $(OBJS)

test: ${TARGET} examples/simple.graphml
	./${TARGET} examples/simple.graphml
	dot -Teps simple_poly.dot -o simple_poly.eps

clean:
	@rm -f $(TARGET)
	@rm -f *.core
	@rm -f *.d *.o
	@rm -f *.dot *.eps
