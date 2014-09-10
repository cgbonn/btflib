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
% Most of Bonn University's BTF binary formats share this common header.
function write_common_header(fid, signature, meta)
    fwrite(fid, signature, 'char');
    
    write_meta_data(fid, meta);
    write_bidir_sampling(fid, meta);
    
    % store spatial extents of textures (in pixels)
    fwrite(fid, meta.width, 'uint32');
    fwrite(fid, meta.height, 'uint32');
    
    % store rotations (if any) of local coordinate frames
    if signature(end) == 'R'
        fwrite(fid, meta.num_rotations, 'uint32');
        if meta.num_rotations > 0
            fwrite(fid, meta.rotations', 'single');
        end
    end
end