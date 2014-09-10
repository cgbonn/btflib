btflib
======

A library providing basic functionality for Bidirectional Texture Functions


dependencies
------------

If you want to use the MATLAB code, you need to download the following external dependencies:

*[Matlab IEEE 754r Half Precision floating point converter](http://www.mathworks.de/matlabcentral/fileexchange/23173-ieee-754r-half-precision-floating-point-converter)*

Please notice that compiling halfprecision.m on 64bit Linux systems requires some changes to the files:

*   in halfprecision.c add the following directives (replacing the existing #defines):

        #include <stdint.h>
        #define INT16_TYPE int16_t
        #define UINT16_TYPE uint16_t
        #define INT32_TYPE int32_t
        #define UINT32_TYPE uint32_t

*   in halfprecision.m replace "mex(cname);" with "eval(['mex -O CFLAGS="\$CFLAGS -std=c99" ', cname]);"
*   you should now be able to compile halfprecision.c by calling halfprecision() without any arguments
