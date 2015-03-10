% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-03-10
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
% Most of Bonn University's BTF binary formats share this common header.
function meta = read_common_header(fid, signature)
    meta = read_meta_data(fid);
    meta = read_bidir_sampling(fid, meta);

    % read spatial extents of textures (in pixels)
    meta.width = fread(fid, 1, 'uint32');
    meta.height = fread(fid, 1, 'uint32');

    % read rotations (if any) of local coordinate frames
    if strcmp(signature(end),'R')
        meta.num_rotations = double(fread(fid, 1, '*uint32'));
        if (meta.num_rotations > 0)
            meta.rotations = read_matrix(fid, '*float', meta.num_rotations, 3, 3);
        end
    else
        meta.num_rotations = 0;
    end
end
