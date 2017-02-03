/*************************************************************************
 * This code is part of Matlab Toolbox.
 * Copyright 2017 University of Bonn
 *
 * authors:
 *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
 *
 * file creation date: 2017-01-25q
 *
 *************************************************************************
 
 Mex file for querying meta data from an OpenEXR image file.
 */
#include <vector>

#include <ImfInputFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfChannelList.h>
#include <ImfPixelType.h>
#include <Iex.h>

#include <mex.h>

#define NUMBER_OF_FIELDS (sizeof(field_names)/sizeof(*field_names))

using namespace Imf;
using namespace Imath;
using namespace Iex;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	if (nrhs != 1 || nlhs > 1 || !mxIsChar(prhs[0]) || mxGetM(prhs[0]) != 1) {
        mexErrMsgTxt("Usage: meta = exr.query(path_to_exr_file)");
    }
    
	char *filename = mxArrayToString(prhs[0]);

	try {
        InputFile file(filename);
        Box2i dw = file.header().dataWindow();

        int height, width, num_channels;
        width = dw.max.x - dw.min.x + 1;
        height = dw.max.y - dw.min.y + 1;
        
        // get channel information from header
        std::vector<std::string> vec_chan_names;
        std::vector<std::string> vec_chan_types;
        const ChannelList &channels = file.header().channels();
        for (ChannelList::ConstIterator i = channels.begin(); i != channels.end(); ++i) {
            vec_chan_names.push_back(i.name());
            
            PixelType type = i.channel().type;
			std::string t = (type == UINT) ? "uint" : ((type == HALF) ? "half" : "float");
            vec_chan_types.push_back(t);
        }
        num_channels = vec_chan_names.size();
        
		const StringAttribute *comments =
            file.header().findTypedAttribute <StringAttribute> ("comments");

        // create meta struct
        mwSize dims[2] = {1, 1};
		const char *field_names[] = {"width", "height", "num_channels", 
            "channel_names", "channel_types", "comments"};
        plhs[0] = mxCreateStructArray(2, dims, NUMBER_OF_FIELDS, field_names);
        
        // set struct fields
        mxArray* field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field_value) = width;
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "width"), field_value);
        
        field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field_value) = height;
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "height"), field_value);
        
        field_value = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(field_value) = num_channels;
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "num_channels"), field_value);
        
        // create cell arrays of strings for the channel names and types
        mxArray* cell_array_ptr_channel_names = mxCreateCellMatrix((mwSize)num_channels, 1);
        mxArray* cell_array_ptr_channel_types = mxCreateCellMatrix((mwSize)num_channels, 1);
        for (int i = 0; i < num_channels; i++) {
            mxSetCell(cell_array_ptr_channel_names, i, mxCreateString(vec_chan_names[i].c_str()));
            mxSetCell(cell_array_ptr_channel_types, i, mxCreateString(vec_chan_types[i].c_str()));
        }
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "channel_names"), cell_array_ptr_channel_names);
        mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "channel_types"), cell_array_ptr_channel_types);
        
        if (comments) {
            mxSetFieldByNumber(plhs[0], 0, mxGetFieldNumber(plhs[0], "comments"), 
                    mxCreateString(comments->value().c_str()));
        }
	} catch (const std::exception &exc) {
		mexErrMsgTxt(exc.what());
	}

	// Free the memory for the string
	mxFree(filename);

	return;
} 
