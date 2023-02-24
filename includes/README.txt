
%compile with
mex trace_flatmap_stack.cc  -largeArrayDims -lgomp CXXFLAGS=" -O3   -Wfatal-errors  -std=c++11 -fopenmp-simd -fopenmp  -fPIC -march=native"
