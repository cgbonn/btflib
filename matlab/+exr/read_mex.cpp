/*************************************************************************
 * This code is part of Matlab Toolbox.
 * Copyright 2017 University of Bonn
 *
 * authors:
 *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
 *
 * file creation date: 2017-01-25
 *
 *************************************************************************
 
 Mex file for reading images in OpenEXR format. Usage:
 
 [image, channel_names] = exr.read_mex(filename, to_single), where
 - the optional argument to_single determines if the pixel values should be
    converted to single precision floats, even if they are stored as unsigned
    integers or as half precision floats in the file
 - image is is a 2D or 3D array of floats or unsigned integers (also for half
    precision floats)
 - channel_names is a cell array of strings holding the names of each channel
 */
#include <cstdint>

#include <ImathBox.h>
#include <ImfChannelList.h>
#include <ImfInputFile.h>
#include <ImfArray.h>
#include <ImfChannelList.h>
#include <ImfPixelType.h>
#include <Iex.h>
#include <half.h>
#include <vector>
#include <fstream>

#include <mex.h>

using namespace Imf;
using namespace Imath;
using namespace Iex;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { 

	if (nrhs < 1 || nrhs > 2 || nlhs < 1 || nlhs > 2) {
		mexErrMsgTxt("Usage: [image, channel_names] = exr.read_mex(filename[, to_single]).");
    }

	if (mxIsChar(prhs[0]) != 1) {
		mexErrMsgTxt("filename argument must be a string.");
    }

	char* filename = mxArrayToString(prhs[0]);
    bool to_single = false;
    if (nrhs > 1) {
        to_single = (bool) mxGetScalar(prhs[1]);
    }

	try {
		InputFile file(filename);
        mxFree(filename);
        Header header = file.header();
        
        // get channel information from header
        std::vector<std::string> channel_names;
        // OpenEXR files can store channels with different data types.
        // Since we want to return one array, we need to determine the smalles
        // and lowest precision data type which stores all channels with the
        // best precision
        PixelType max_pixel_type = UINT;
        const ChannelList &channels = header.channels();
        for (ChannelList::ConstIterator i = channels.begin(); i != channels.end(); i++) {
            channel_names.push_back(i.name());
            
            // the channel with the highest precision data type determines the output precision
            // the ordering of the PixelType enum allows to use max here
            max_pixel_type = std::max(max_pixel_type, i.channel().type);
        }
        
        // output format (HALF & UINT both are represented as uint16 in Matlab)
        mxClassID output_type = max_pixel_type == FLOAT ? mxSINGLE_CLASS : 
            (max_pixel_type == HALF ? mxUINT16_CLASS : mxUINT32_CLASS);
        
        // enforce conversion to single precision floats?
        if (to_single) {
            max_pixel_type = FLOAT;
            output_type = mxSINGLE_CLASS;
        }

        // set up output dimensions (warning: this ignores all datawindow related stuff)
        Box2i dw = header.dataWindow();
		int width = dw.max.x - dw.min.x + 1;
		int height = dw.max.y - dw.min.y + 1;
        unsigned num_channels = (unsigned)channel_names.size();
        
        int dims[3];
		dims[0] = height; 
		dims[1] = width; 
		dims[2] = num_channels;
        
        // reserve output array
        plhs[0] = mxCreateNumericArray(3, dims, output_type, mxREAL);
        
        FrameBuffer fb;
        
        // set up slices per channel with the determined output format
        if (max_pixel_type == UINT) {
            uint32_t* img = (uint32_t*) mxGetData(plhs[0]);
            for(unsigned int i = 0; i < num_channels; i++) {
                fb.insert(channel_names[i].c_str(),
                        Slice(UINT, (char*) &img[i * height * width],
                        sizeof(uint32_t) * height, sizeof(uint32_t)));
            }
        } else if (max_pixel_type == HALF) {
            half* img = (half*) mxGetData(plhs[0]);
            for(unsigned int i = 0; i < num_channels; i++) {
                fb.insert(channel_names[i].c_str(),
                        Slice(HALF, (char*) &img[i * height * width],
                        sizeof(half) * height, sizeof(half)));
            }
        } else {
            float* img = (float*) mxGetData(plhs[0]);
            for(unsigned int i = 0; i < num_channels; i++) {
                fb.insert(channel_names[i].c_str(),
                        Slice(FLOAT, (char*) &img[i * height * width],
                        sizeof(float) * height, sizeof(float)));
            }
        }
        
        // perform the actual copy process
        file.setFrameBuffer(fb);
        file.readPixels(dw.min.y, dw.max.y);
        
        // channel names
		if (nlhs > 1) {
            dims[0] = 1;
            dims[1] = num_channels;
            
            plhs[1] = mxCreateCellArray(2, dims);
            for (unsigned int i = 0; i < num_channels; i++) {
                mxSetCell(plhs[1], i, mxDuplicateArray(mxCreateString(channel_names[i].c_str())));
            }
		}
	} catch (const std::exception &exc) {
		mxFree(filename);
		mexErrMsgTxt(exc.what());
	}

	return;
} 
