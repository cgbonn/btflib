btflib
======

A library providing basic functionality for Bidirectional Texture Functions


dependencies
------------

If you want to use the MATLAB code, you need to compile the bundled IEEE 754r Half Precision floating point library:

*   change to the subfolder *matlab/thirdparty/half_precision/*
*   from the MATLAB command window call halfprecision() without any arguments
*   this should automatically call the mex compiler and if no errors are produced, you can test the library with
*   halfprecision(halfprecision(single(0.1)), 'single')
*   this should output 0.1 again
*   you can now add the subfolder *matlab/thirdparty/half_precision/* to your MATLAB search path: addpath('matlab/thirdparty/half_precision/') (from this repository's root folder)

The following changes have been applied to the original code from http://www.mathworks.de/matlabcentral/fileexchange/23173-ieee-754r-half-precision-floating-point-converter to make it compile on 64bit Linux systems:

*   in halfprecision.c add the following directives (replacing the existing #defines):

        #include <stdint.h>
        #define INT16_TYPE int16_t
        #define UINT16_TYPE uint16_t
        #define INT32_TYPE int32_t
        #define UINT32_TYPE uint32_t

*   in halfprecision.m replace "mex(cname);" with "eval(['mex -O CFLAGS="\$CFLAGS -std=c99" ', cname]);"
