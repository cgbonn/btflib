% *************************************************************************
% * Copyright 2017 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * file creation date: 2017-01-25
% *
% * This file is part of btflib.
% *
% * btflib is free software: you can redistribute it and/or modify it under
% * the terms of the GNU Lesser General Public License as published by the
% * Free Software Foundation, either version 3 of the License, or (at your
% * option) any later version.
% *
% * btflib is distributed in the hope that it will be useful, but WITHOUT
% * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
% * License for more details.
% *
% * You should have received a copy of the GNU Lesser General Public
% * License along with btflib.  If not, see <http://www.gnu.org/licenses/>.
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
    
    mpath = fileparts(mfilename('fullpath'));
    structmex = dir(fullfile(mpath, ['read_mex.', mexext]));
    structcpp = dir(fullfile(mpath, 'read_mex.cpp'));
    if isempty(structmex) || structmex.datenum < structcpp.datenum
        % compile mex file if it cannot be found
        warning('exrread:mex_outdated', ...
            ['the mex file exr.read_mex is outdated ', ...
            'or non existant and needs to compiled.']);
        if isunix
            mex(fullfile(mpath, 'read_mex.cpp'), '-I/usr/include/OpenEXR', ...
                '-lIlmImf', '-lHalf', '-outdir', mpath);
        else
            mex(fullfile(mpath, 'read_mex.cpp'), ...
                '-Lc:\path\to\OpenEXR\lib\', ...
                '-Ic:\path\to\OpenEXR\include', ...
                '-Ic:\path\to\IlmBase\', ...
                '-Lc:\path\to\zlib\lib\', ...
                '-lIlmImf', '-lIex', '-lImath', '-lHalf', '-lzlib', ...
                '-outdir', mpath);
        end
    end
    
    [im, channel_names] = exr.read_mex(filename, to_single);
    
    if numel(channel_names) == 3 && all(cellfun(@isequal, channel_names, {'B', 'G', 'R'}))
        im = im(:, :, 3 : -1 : 1);
        channel_names = channel_names(end : -1 : 1);
    end
end
