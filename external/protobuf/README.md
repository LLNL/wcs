## Building the protocol buffers library
   Here, we describe two ways to build this library.
 + One is to leave the job to the WCS build system. When WCS cannot find an
   existing installation of protobuf library, it will automatically download
   the source of the library, build and install it.
   The whole process of building WCS consists of two steps in this case. The
   first step will install the protobuf. The second step will build WCS using
   the protobuf installed.
 + To build a stand-alone copy of the libarary under a directory out of
   the WCS source tree, users can use the `CMakeLists.txt` provided here
   with the option `-DCMAKE_INSTALL_PREFIX=<installation-path>`
   No need to run `make install`. It is done with `make`.
