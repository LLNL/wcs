BASE_DIR = .
CUR_DIR := $(BASE_DIR)/utils
include $(CUR_DIR)/../config.mk
CPPFLAGS += -I$(BASE_DIR)
LIBS += -lboost_filesystem -lboost_system

.PHONY: all target_dep target_rule test clean

CUR_TARGET := $(CUR_DIR)/graphml2dot
SRCS := $(CUR_DIR)/graph_factory.cpp $(CUR_DIR)/exception.cpp $(CUR_DIR)/trace.cpp $(CUR_DIR)/trace_ssa.cpp $(CUR_TARGET).cpp

OBJS := $(addsuffix .o, $(basename $(SRCS)))

WCM_OBJS := $(BASE_DIR)/reaction_network/*.o $(CUR_DIR)/graph_factory.o $(CUR_DIR)/exception.o $(CUR_DIR)/trace.o $(CUR_DIR)/trace_ssa.o


${CUR_TARGET}: $(SRCS) $(WCM_OBJS)
	@echo "building $(CUR_TARGET)"
	${CXX} ${CXXFLAGS} ${CPPFLAGS} -o $@ $@.cpp $(WCM_OBJS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

%.d: %.cpp
	@echo "$@ `$(CXXDEP) $(CXXFLAGS) $(CPPFLAGS) -MM $<`" > $@

target_dep: $(addsuffix .d, $(basename $(SRCS)))

target_rule: $(OBJS)

ifneq ($(MAKECMDGOALS),clean)
    -include $(OBJS:.o=.d)
endif

test: $(CUR_TARGET) $(BASE_DIR)/tests/example/simple.graphml
	$(CUR_TARGET) $(BASE_DIR)/tests/example/simple.graphml
	dot -Teps simple.dot -o simple.eps

clean:
	@rm -f $(CUR_TARGET)
	@rm -f $(CUR_DIR)/*.core
	@rm -f $(CUR_DIR)/*.d $(CUR_DIR)/*.o
	@rm -f $(CUR_DIR)/*.dot $(CUR_DIR)/*.eps
