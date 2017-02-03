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
% Query meta data from an OpenEXR image file.
function meta = query(filename)
    if exist('+exr/query_mex', 'file') ~= 3
        % compile mex file if it cannot be found
        mpath = fileparts(mfilename('fullpath'));
        if isunix
            mex(fullfile(mpath, 'query_mex.cpp'), '-I/usr/include/OpenEXR', ...
                '-lIlmImf', '-lHalf', '-outdir', mpath);
        else
            mex(fullfile(mpath, 'query_mex.cpp'), ...
                '-Lc:\path\to\OpenEXR\lib\', ...
                '-Ic:\path\to\OpenEXR\include', ...
                '-Ic:\path\to\IlmBase\', ...
                '-Lc:\path\to\zlib\lib\', ...
                '-lIlmImf', '-lIex', '-lImath', '-lHalf', '-lzlib');
        end
    end
    
    meta = exr.query_mex(filename);
end
