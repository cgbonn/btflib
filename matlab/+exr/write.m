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
% Function for writing images in OpenEXR format. Usage:
%
% exr.write(image, filename[, write_half[, channel_names[, compression]]]),
% where
% - image is is a 2D or 3D array of floats, uint32s or uint16s (for half
%   precision floats) 
% - the optional argument precision enforces the data to be interpreted as the
%   specified data type, possible values are 'single', 'half' or 'uint'
% - channel_names is a cell array of strings holding the names of each channel
% - compression is a string with one of the following values:
%   - none:  no compression
%   - rle:   run length encoding
%   - zips:  zlib compression, one scan line at a time
% 	- zip:   zlib compression, in blocks of 16 scan lines
% 	- piz:   piz-based wavelet compression
% 	- pxr24: lossy 24-bit float compression
% 	- b44:   lossy 4-by-4 pixel block compression, fixed compression rate
% 	- b44a:  lossy 4-by-4 pixel block compression, flat fields are compressed
% 	         more
function write(im, filename, precision, channel_names, compression)
    if ~exist('precision', 'var') || isempty(precision)
        if isa(im, 'uint16')
            precision = 'half';
        elseif isa(im, 'uint32')
            precision = 'uint';
        else
            precision = 'single';
        end
    end
    
    if ~ismember(precision, {'uint', 'uint32', 'half', 'single', 'float'})
        error('exrwrite:precision', ...
            'The precision argument must be one of ''uint'', ''half'' or ''single''.');
    end
    
    % scalars are so much nicer than strings :)
    if strcmpi(precision, 'uint') || strcmpi(precision, 'uint32')
        precision = 0;
    elseif strcmpi(precision, 'half')
        precision = 1;
    else
        precision = 2;
    end
    
    % OpenEXR only supports single floats as highest precision data type
    if isa(im, 'double')
        im = single(im);
    end
    
    if exist('channel_names', 'var') && ~iscell(channel_names)
        if ischar(channel_names) && numel(channel_names) == size(im, 3)
            channel_names = num2cell(channel_names);
        else
            channel_names = {channel_names};
        end
    end
    
    if ~exist('channel_names', 'var')
        if size(im, 3) == 1
            channel_names = {'L'};
        elseif size(im, 3) == 3
            channel_names = {'R', 'G', 'B'};
        else
            channel_names = utils.sprintf2('C%d', 1 : size(im, 3));
        end
    end
    
    % ensure image is handed to c++ in the right format
    if precision == 0 && ~isa(im, 'uint32')
        im = uint32(im);
    elseif precision == 1 && ~isa(im, 'uint16')
        if isfloat(im)
            im = halfprecision(im);
        elseif isa(im, 'int16')
            warning('exrwrite:halfprecision', ...
                'converting int16 to uint16');
            im = typecast(im, 'uint16');
        else
            error('exrwrite:nohalfprecision', ...
                'image is not in float or uint16 format');
        end
    elseif precision == 2 && ~isa(im, 'single')
        im = single(im);
    end
    
    if ~exist('compression', 'var') || isempty(compression)
        compression = 'none';
    end
    
    if ~(ismember(lower(compression), ...
            {'no', 'none', 'rle', 'zips', 'zip', 'piz', 'pxr24', 'b44', 'b44a'}) ...
            || isscalar(compression) && ismember(compression, 0 : 7))
        error('exrwrite:compression_format', ...
            ['unknown compression format (%s), please select one of: ', ...
            'none, rle, zips, zip, piz, pxr24, b44, b44a'], compression);
    end
    
    if ischar(compression)
        switch lower(compression)
            case {'no', 'none'}
                compression = 0;
            case 'rle'
                compression = 1;
            case 'zips'
                compression = 2;
            case 'zip'
                compression = 3;
            case 'piz'
                compression = 4;
            case 'pxr24'
                compression = 5;
            case 'b44'
                compression = 6;
            case 'b44a'
                compression = 7;
            otherwise
                error('exrwrite:unknown_compression', 'unknown compression method');
        end
    end
    
    % check if compiled mex file exists or is outdated
    mpath = fileparts(mfilename('fullpath'));
    structmex = dir(fullfile(mpath, ['write_mex.', mexext]));
    structcpp = dir(fullfile(mpath, 'write_mex.cpp'));
    if isempty(structmex) || structmex.datenum < structcpp.datenum
        % compile mex file if it cannot be found or if it is outdated
        warning('exrquery:mex_outdated', ...
            ['the mex file exr.write_mex is outdated ', ...
            'or non existant and needs to compiled.']);
        if isunix
            mex(fullfile(mpath, 'write_mex.cpp'), '-I/usr/include/OpenEXR', ...
                '-lIlmImf', '-lHalf', '-outdir', mpath);
        else
            mex(fullfile(mpath, 'write_mex.cpp'), ...
                '-Lc:\path\to\OpenEXR\lib\', ...
                '-Ic:\path\to\OpenEXR\include', ...
                '-Ic:\path\to\IlmBase\', ...
                '-Lc:\path\to\zlib\lib\', ...
                '-lIlmImf', '-lIex', '-lImath', '-lHalf', '-lzlib', ...
                '-outdir', mpath);
        end
    end
    
    exr.write_mex(im, filename, precision(1), channel_names, compression);
end
