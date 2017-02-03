% *************************************************************************
% * This code is part of Matlab Toolbox.
% * Copyright 2017 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * file creation date: 2017-01-25
% *
% *************************************************************************
%
% Function for reading images in OpenEXR format. Usage:
%
% [image, channel_names] = exr.read(filename[, to_single]), where
% - the optional argument to_single determines if the pixel values should be
%   converted to single precision floats, even if they are stored as unsigned
%   integers or as half precision floats in the file
% - image is is a 2D or 3D array of floats or unsigned integers (also for half
%   precision floats)
% - channel_names is a cell array of strings holding the names of each channel
function [im, channel_names] = read(filename, to_single)
    if ~exist('to_single', 'var')
        to_single = true;
    end
    
    if exist('+exr/read_mex', 'file') ~= 3
        % compile mex file if it cannot be found
        mpath = fileparts(mfilename('fullpath'));
        if isunix
            mex(fullfile(mpath, 'read_mex.cpp'), '-I/usr/include/OpenEXR', ...
                '-lIlmImf', '-lHalf', '-outdir', mpath);
        else
            mex(fullfile(mpath, 'read_mex.cpp'), ...
                '-Lc:\path\to\OpenEXR\lib\', ...
                '-Ic:\path\to\OpenEXR\include', ...
                '-Ic:\path\to\IlmBase\', ...
                '-Lc:\path\to\zlib\lib\', ...
                '-lIlmImf', '-lIex', '-lImath', '-lHalf', '-lzlib');
        end
    end
    
    [im, channel_names] = exr.read_mex(filename, to_single);
    
    if numel(channel_names) == 3 && all(cellfun(@isequal, channel_names, {'B', 'G', 'R'}))
        im = im(:, :, 3 : -1 : 1);
        channel_names = channel_names(end : -1 : 1);
    end
end
