% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-12-22
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
% Given one or several vectors and optionally a single or the same number
% of normal vectors, this function computes the reflection direction(s),
% i.e. the vector with inclination angle rotated by 180 degrees around the
% normal vector.
function reflected_dirs = reflect(dirs, normals)
    if ~exist('normals', 'var')
        normals = [0; 0; 1];
    end
    
    % bring input to 3 x N format
    s_in = size(dirs);
    dirs = arrange(dirs);
    normals = arrange(normals);
    
    nd = size(dirs, 2);
    nn = size(normals, 2);
    
    % repeat normal if necessary
    if nn == 1
        normals = repmat(normals, 1, nd);
    elseif nd ~= nn
        error('Please input the same number of directions and normals, or only a single normal.');
    end
    
    % compute the reflection
    reflected_dirs = repmat(2 .* dot(normals, dirs, 1), [3, 1]) .* normals - dirs;
    
    % bring output into the same shape as the input
    if numel(s_in) == 2 && s_in(2) == 3 && nd ~= 3
        reflected_dirs = reflected_dirs';
    else
        d = find(s_in == 3, 1);
        reflected_dirs = reshape(reflected_dirs, [3, s_in(setdiff(1 : numel(s_in), d))]);
        reflected_dirs = ipermute(reflected_dirs, [d, setdiff(1 : numel(s_in), d)]);
    end
    
    function mat = arrange(mat)
        s = size(mat);
        if s(1) ~= 3
            dim = find(s == 3, 1);
            mat = permute(mat, [dim, setdiff(1 : numel(s), dim)]);
        end
        mat = reshape(mat, 3, []);
    end
end
