### A note for setting up performance profiling

 Often, we are not interested in the cost of a certain region of the code
 but only in the performance of a specific region that we target to optimize.
 In such a case, we need to indicate the code region to profile such that a
 profiling tool can avoid measuring uninterested regions of the code, and
 focus on collecting data from where we are interested in.

 + [**VTune**]( https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/vtune-profiler.html)
    To mark the code region to profile, we compile the code with a header
    and link with the static library of the tool such that we can insert the
    calls to start and stop profiling, which are
    ` __itt_resume()` and `__itt_pause()`.

    On LC, load the vtune module, and then find out the location of the tool by
      `module load vtune; module show vtune`

    Use the following cmake options when building WCS:
      `-DWCS_WITH_VTUNE:BOOL=ON -DVTUNE_DIR:PATH=...`

 + [**HPCToolkit**]( http://hpctoolkit.org)
    To mark the code region to profile, we compile the code with the header
    and link with the dynamic library of the tool such that we can insert the 
    calls to start and stop profiling, which are
    ` hpctoolkit_sampling_start()` and ` hpctoolkit_sampling_stop()`.
    We use the dynamic one instead of static one to avoid using the linker 
    wrapper provided by the tool for instrumentation.

    On LC, load the hpctoolkit module, and find out the location of the tool by
      `module load hpctoolkit; module show hpctoolkit`

    Use the following cmake options when building WCS:
      `-DWCS_WITH_HPCTOOLKIT:BOOL=ON -DHPCTOOLKIT_DIR:PATH=...`

 + [**gprof**]( https://hpc.llnl.gov/software/development-environment-software/gprof)
    To compile the code for generating the measurement file that can be
    analyzed by gprof, use the cmake option `-DWCS_GPROF:BOOL=ON`.
    It will add `-pg` compiler option.
    We are not able to define a code region to profile using gprof. However,
    we can selectively disable an uninteresting region of the code as described
    below using the macro definition `WCS_PERF_PROF`

 + Misc.
    See `src/reaction.cpp` to see an example of how the start and stop is
    instrumented. Also, see the generated `wcs_config.hpp` that includes
    the headers for the profiling tools.

    If any of the cmake options `WCS_PERF_PROF`, `WCS_WITH_HPCTOOLKIT`
    or `WCS_WITH_VTUNE` is used, the macro definition `WCS_PERF_PROF`
    can be used to selectively disable a region of code to avoid profiling that
    part of the code.

    ```
    #ifndef WCS_PERF_PROF
     lines_to_be_disabled_during_profiling ....
    #endif
    ```
