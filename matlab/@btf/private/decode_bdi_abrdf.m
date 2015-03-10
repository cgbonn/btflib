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
% Decode a single ABRDF from a bidirectional image file (BDI).
function abrdf = decode_bdi_abrdf(obj, x, y)
    [chunk_index, chunk_abrdf_index] = obj.get_bdi_chunk_index(x, y);
    
    if chunk_index == -1
        abrdf = zeros(obj.meta.nL, obj.meta.nV, obj.meta.num_channels, 'single');
        return;
    end
    
    if obj.data.chunks_buffered(chunk_index)
        abrdf = obj.data.chunks((chunk_abrdf_index - 1) * obj.meta.abrdf_size + 1: ...
            chunk_abrdf_index * obj.meta.abrdf_size, chunk_index);
    else
        obj.data.fid = fopen(obj.meta.file_name, 'r');
        abrdf = obj.get_bdi_chunk([], x, y);
        fclose(obj.data.fid);
    end
    
    % convert half precision floats to single precision floats & rearrange data
    abrdf = permute(reshape(halfprecision(abrdf, 'single'), ...
        obj.meta.num_channels, obj.meta.nL, obj.meta.nV), [2, 3, 1]);
end
