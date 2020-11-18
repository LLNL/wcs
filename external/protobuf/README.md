## Building the protocol buffers library
   Here, we describe two ways to build this library.
 + One is to leave the job to the WCS build system. By default, WCS
   will automatically download the source of the library. Then, it will be
   built and installed along with the rest of WCS project.
 + To build a stand-alone copy of the libarary under a directory out of
   the WCS source tree, users can use the `CMakeLists.txt` provided here
   with the option `-DCMAKE_INSTALL_PREFIX=<installation-path>`
