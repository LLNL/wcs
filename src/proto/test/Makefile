PROTO_DIR = ../../../install

PROTO_INC_DIR = $(PROTO_DIR)/include
PROTO_LIB_DIR = $(PROTO_DIR)/lib64
PROTOC = $(PROTO_DIR)/bin/protoc
CXXFLAGS = -Wall -std=c++17 -m64 -I../.. -I$(PROTO_INC_DIR)

ifndef PROTO_DIR
  $(error Set PROTO_DIR variable to the path of the protocol buffer installation)
endif

# In case of using clang, pick the right std library that is compatible with
# the one used in building the protocol buffer

ifeq ($(CXX), clang++)
    #STDLIB = --stdlib=libc++
    STDLIB = --gcc-toolchain=/usr/tce/packages/gcc/gcc-8.1.0
endif

LDFLAGS = $(STDLIB) -m64 -Wl,-rpath=$(PROTO_LIB_DIR),--no-undefined -L$(PROTO_LIB_DIR) -lprotobuf

all: $(PROTOC) wcs_proto_test wcs_params.protobin

$(PROTOC):
	@echo "Checking if '$(PROTOC)' exists ..."
	@test -x $(PROTOC) || { echo "'$(PROTOC)' does not exist!"; echo "Set PROTO_DIR to a correct path."; }
	@exit 1

wcs_params.pb.cc: $(PROTOC) ../wcs_params.proto
	$(PROTOC) --proto_path=.. --cpp_out=. wcs_params.proto

wcs_params.pb.o: wcs_params.pb.cc
	$(CXX) $(CXXFLAGS) $(STDLIB) -c wcs_params.pb.cc -o wcs_params.pb.o

wcs_proto_test: wcs_proto_test.cpp wcs_params.pb.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o wcs_proto_test wcs_proto_test.cpp wcs_params.pb.o

# Convert the prototext file into a binary file
wcs_params.protobin: $(PROTOC) wcs_params.prototext
	$(PROTOC) --proto_path=.. --encode="wcs_proto.WCS_Params" wcs_params.proto < wcs_params.prototext > wcs_params.protobin

# Try reading both the prototext input and its binary counterpart
test: wcs_proto_test wcs_params.protobin
	./wcs_proto_test wcs_params.prototext > readtext.txt
	./wcs_proto_test wcs_params.protobin 1 > readbin.txt
	@diff readtext.txt wcs_params.prototext
	@diff readbin.txt wcs_params.prototext
	
clean:
	@rm -f *.o wcs_proto_test
	@rm -f wcs_params.protobin
	@rm -f wcs_params.pb.cc  wcs_params.pb.h
	@rm -f readtext.txt readbin.txt
