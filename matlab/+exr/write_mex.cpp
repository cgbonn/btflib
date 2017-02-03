/**************************************************************************
 * This code is part of Matlab Toolbox.
 * Copyright 2017 University of Bonn
 *
 * authors:
 *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
 *
 * file creation date: 2017-01-25
 *
 *************************************************************************

 Mex file for writing images in OpenEXR format. Usage:

 exr.write_mex(image, filename[, write_half[, channel_names]]), where
 - image is is a 2D or 3D array of floats, uint32s or uint16s (for half
   precision floats) 
 - the optional argument precision enforces the data to be interpreted as the
   specified data type, possible values are 'single', 'half' or 'uint'
 - channel_names is a cell array of strings holding the names of each channel
 */
#include <algorithm>
#include <cstdint>
#include <vector>

#include <half.h>
#include <ImathBox.h>
#include <ImfHeader.h>
#include <ImfOutputFile.h>
#include <ImfArray.h>
#include <ImfChannelList.h>
#include <ImfPixelType.h>
#include <Iex.h>

#include <mex.h>
#include <matrix.h>

using namespace Imf;
using namespace Imath;
using namespace Iex;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // check & parse inputs
    if (nlhs != 0) {
		mexErrMsgTxt("Function does not return any outputs.");
    }

	if (nrhs != 4) {
		mexErrMsgTxt("Usage: exr.write_mex(image, filename, precision, channel_names);");
    }
    
    if (!mxIsChar(prhs[1])) {
		mexErrMsgTxt("Second input argument must be a string.");
    }
    
    if (!mxIsCell(prhs[3]) || mxIsEmpty(prhs[3])) {
        mexErrMsgTxt("input must either be M x N x 3 or M x N x P and a cell array of P strings specifying the channel names.");
    }
    
    char *filename = mxArrayToString(prhs[1]);
    int precision = (int) mxGetScalar(prhs[2]); // 0: UINT, 1: HALF, 2: FLOAT
    const mxArray* mx_channels = prhs[3];
    mwSize num_channel_names = mxGetNumberOfElements(mx_channels);

	int ndims = mxGetNumberOfDimensions(prhs[0]);
	if (ndims < 2 || ndims > 3) {
		mexErrMsgTxt("First input argument must be a 2D or 3D array.");
    }
    const mwSize* dims = mxGetDimensions(prhs[0]);
    int height = dims[0];
    int width = dims[1];
    int num_channels = 1;
    if (ndims == 3 && mxIsCell(prhs[3])) {
        num_channels = dims[2];
    }
    
    if (num_channels != num_channel_names) {
        mexErrMsgTxt("Number of image channels must match number of channel names!");
    }
    
    if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS 
            && mxGetClassID(prhs[0]) != mxUINT16_CLASS 
            && mxGetClassID(prhs[0]) != mxUINT32_CLASS) {
		mexErrMsgTxt("First input argument must be in single or uint16 precision.");
    }
    
    // conversion from lower to higher precision doesn't make sense (except for uint16 -> uint32)
    if (precision == 2 && (mxGetClassID(prhs[0]) == mxUINT16_CLASS 
            || mxGetClassID(prhs[0]) == mxUINT32_CLASS)) {
        mexErrMsgTxt("If the image array is in uint16 (or half) or uint32 format, precision must be set to 'half' or 'uint'.");
    }

	try {
        float* img = (float*) mxGetData(prhs[0]);
        half* img_half = (half*) mxGetData(prhs[0]);
        unsigned int* img_uint = (unsigned int*) mxGetData(prhs[0]);
        int do_cleanup = -1;
        
        // convert to uint / half precision if requested but image is in different format
        if (precision == 0 && mxGetClassID(prhs[0]) != mxUINT32_CLASS) {
            img_uint = new unsigned int[height * width * num_channels];
            do_cleanup = 0;
            for (int i = 0; i < height * width * num_channels; i++) {
                img_uint[i] = (unsigned int) img[i];
            }
        } else if (precision == 1 && mxGetClassID(prhs[0]) != mxUINT16_CLASS) {
            img_half = new half[height * width * num_channels];
            do_cleanup = 1;
            for (int i = 0; i < height * width * num_channels; i++) {
                img_half[i] = (half) img[i];
            }
        }
        
        Header header(width, height);

        std::vector<std::string> channel_names(num_channels);
        for (int i = 0; i < num_channels; i++) {
            const mxArray* cell_ptr = mxGetCell(mx_channels, i);
            mwSize buflen = mxGetN(cell_ptr) * sizeof(mxChar) + 1;
            char* c_array = (char*) mxMalloc(buflen);
            int status = mxGetString(cell_ptr, c_array, buflen);
            if (status == 0) {
                channel_names[i] = std::string(c_array);
                header.channels().insert(c_array, 
                        Channel(precision == 0 ? UINT : (precision == 1 ? HALF : FLOAT)));
            } else {
                mexErrMsgTxt("Channel names must be specified as cell array of strings!");
            }
            mxFree(c_array);
        }

        OutputFile file(filename, header);
        FrameBuffer framebuffer;
        for (int i = 0; i < num_channels; i++) {
            char* p_slice = (char*) &img[i * height * width];
            int x_stride = sizeof(float) * height;
            int y_stride = sizeof(float);
            if (precision == 0) {
                p_slice = (char*) &img_uint[i * height * width];
                x_stride = sizeof(unsigned int) * height;
                y_stride = sizeof(unsigned int);
            } else if (precision == 1) {
                p_slice = (char*) &img_half[i * height * width];
                x_stride = sizeof(half) * height;
                y_stride = sizeof(half);
            }
            framebuffer.insert(channel_names[i], 
                    Slice(precision == 0 ? UINT : (precision == 1 ? HALF : FLOAT),
                    p_slice, x_stride, y_stride));
        }

        file.setFrameBuffer(framebuffer);
        file.writePixels(height);
        
        mxFree(filename);
        if (do_cleanup == 0) {
            delete[] img_uint;
            img_uint = NULL;
        } else if (do_cleanup == 1) {
            delete[] img_half;
            img_half = NULL;
        }
	} catch (const std::exception &exc) {
		mxFree(filename);
		mexErrMsgTxt(exc.what());
	}

	return;
}
