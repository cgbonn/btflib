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
% Function for writing images in OpenEXR format. Usage:
%
% exr.write(image, filename[, write_half[, channel_names]]), where
% - image is is a 2D or 3D array of floats, uint32s or uint16s (for half
%   precision floats) 
% - the optional argument precision enforces the data to be interpreted as the
%   specified data type, possible values are 'single', 'half' or 'uint'
% - channel_names is a cell array of strings holding the names of each channel
function write(im, filename, precision, channel_names)
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
        error('The precision argument must be one of ''uint'', ''half'' or ''single''.');
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
    
    if ~exist('channel_names', 'var')
        if size(im, 3) == 1
            channel_names = {'L'};
        elseif size(im, 3) == 3
            channel_names = {'R', 'G', 'B'};
        else
            channel_names = utils.sprintf2('C%d', 1 : size(im, 3));
        end
    end
    
    if exist('+exr/write_mex', 'file') ~= 3
        % compile mex file if it cannot be found
        mpath = fileparts(mfilename('fullpath'));
        if isunix
            mex(fullfile(mpath, 'write_mex.cpp'), '-I/usr/include/OpenEXR', ...
                '-lIlmImf', '-lHalf', '-outdir', mpath);
        else
            mex(fullfile(mpath, 'write_mex.cpp'), ...
                '-Lc:\path\to\OpenEXR\lib\', ...
                '-Ic:\path\to\OpenEXR\include', ...
                '-Ic:\path\to\IlmBase\', ...
                '-Lc:\path\to\zlib\lib\', ...
                '-lIlmImf', '-lIex', '-lImath', '-lHalf', '-lzlib');
        end
    end
    
    exr.write_mex(im, filename, precision(1), channel_names);
end
