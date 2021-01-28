BOOST_DIR = /usr/tce/packages/boost/boost-1.69.0-mvapich2-2.3-gcc-8.1.0
PROTO_DIR = ../../../install

BOOST_INC_DIR = $(BOOST_DIR)/include
BOOST_LIB_DIR = $(BOOST_DIR)/lib

PROTO_INC_DIR = $(PROTO_DIR)/include
PROTO_LIB_DIR = $(PROTO_DIR)/lib64
PROTOC = $(PROTO_DIR)/bin/protoc

CXXFLAGS = -Wall -std=c++17 -m64 -I../.. -I$(BOOST_INC_DIR) -I$(PROTO_INC_DIR)

ifndef PROTO_DIR
  $(error Set PROTO_DIR variable to the path of the protocol buffer installation)
endif

# In case of using clang, pick the right std library that is compatible with
# the one used in building the protocol buffer

ifeq ($(CXX), clang++)
    #STDLIB = --stdlib=libc++
    STDLIB = --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.1.0
endif

LDFLAGS = $(STDLIB) -m64 -Wl,-rpath=$(BOOST_LIB_DIR),--no-undefined -L$(BOOST_LIB_DIR) -lboost_program_options -Wl,-rpath=$(PROTO_LIB_DIR),--no-undefined -L$(PROTO_LIB_DIR) -lprotobuf

all: $(PROTOC) wcs_proto_progopt_test

$(PROTOC):
	@echo "Checking if '$(PROTOC)' exists ..."
	@test -x $(PROTOC) || { echo "'$(PROTOC)' does not exist!"; echo "Set PROTO_DIR to a correct path."; }
	@exit 1

wcs_params.pb.cc: $(PROTOC) ../wcs_params.proto
	$(PROTOC) --proto_path=.. --cpp_out=. wcs_params.proto

wcs_params.pb.o: wcs_params.pb.cc
	$(CXX) $(CXXFLAGS) $(STDLIB) -c wcs_params.pb.cc -o wcs_params.pb.o

wcs_proto_progopt_test: wcs_proto_progopt_test.cpp wcs_params.pb.o wcs_params.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o wcs_proto_progopt_test wcs_proto_progopt_test.cpp wcs_params.pb.o wcs_params.o

wcs_params.o: ../../wcs_params.cpp ../../wcs_params.hpp
	$(CXX) $(CXXFLAGS) -c ../../wcs_params.cpp

test: wcs_proto_progopt_test
	@echo ----------------------------------------------------
	./wcs_proto_progopt_test --setup-all ../test/wcs_params.prototext --input-model model
	@echo ----------------------------------------------------
	./wcs_proto_progopt_test --setup-all ../test/wcs_params.prototext --input-model model --setup-sim sim_setup.prototext
	@echo ----------------------------------------------------
	./wcs_proto_progopt_test --input-model model --setup-sim sim_setup.prototext
	@echo ----------------------------------------------------
	./wcs_proto_progopt_test --setup-sim sim_setup.prototext --setup-part part_setup.prototext --setup-des des_setup.prototext model
	@echo ----------------------------------------------------

clean:
	@rm -f *.o wcs_proto_progopt_test
	@rm -f wcs_params.pb.cc  wcs_params.pb.h
