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
% Query meta data from an OpenEXR image file.
function meta = query(filename)
    mpath = fileparts(mfilename('fullpath'));
    structmex = dir(fullfile(mpath, ['query_mex.', mexext]));
    structcpp = dir(fullfile(mpath, 'query_mex.cpp'));
    if isempty(structmex) || structmex.datenum < structcpp.datenum
        % compile mex file if it cannot be found
        warning('exrquery:mex_outdated', ...
            ['the mex file exr.query_mex is outdated ', ...
            'or non existant and needs to compiled.']);
        if isunix
            mex(fullfile(mpath, 'query_mex.cpp'), '-I/usr/include/OpenEXR', ...
                '-lIlmImf', '-lHalf', '-outdir', mpath);
        else
            mex(fullfile(mpath, 'query_mex.cpp'), ...
                '-Lc:\path\to\OpenEXR\lib\', ...
                '-Ic:\path\to\OpenEXR\include', ...
                '-Ic:\path\to\IlmBase\', ...
                '-Lc:\path\to\zlib\lib\', ...
                '-lIlmImf', '-lIex', '-lImath', '-lHalf', '-lzlib', ...
                '-outdir', mpath);
        end
    end
    
    meta = exr.query_mex(filename);
end
