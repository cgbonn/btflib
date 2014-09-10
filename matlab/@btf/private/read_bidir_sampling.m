% *************************************************************************
% * Copyright 2014 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2014-09-10
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
% Read bidirectional sampling to a file in one of Bonn University's binary
% formats.
function meta = read_bidir_sampling(fid, meta)
    meta.num_views = fread(fid, 1, '*uint32');
    V = zeros(meta.num_views, 2);
    L = [];
    for vi = 1 : meta.num_views
        V(vi, :) = fread(fid, 2, 'single');
        num_lights = fread(fid, 1, '*uint32');
        if (isempty(L))
            L = fread_matrix(fid, 'single', num_lights, 2);
            meta.num_lights = num_lights;
        elseif (num_lights ~= meta.num_lights)
            error('currently only fixed light hemisphere allowed!!');
        else
            fread(fid,  2 * num_lights, 'single'); %skip light directions
        end
    end
    clear v num_lights;

    meta.L = sph2cart2(L);
    meta.V = sph2cart2(V);
    
    meta.nL = uint32(size(L, 1));
    meta.nV = uint32(size(V, 1));
end