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
% Given two texture coordinates, determine the chunk index and the offset inside
% the chunk for the ABRDF corresponding to that pixel.
function [chunk_index, chunk_abrdf_index] = get_bdi_chunk_index(obj, x, y)
    logical_abrdf_index = (y - 1) * obj.meta.width + x;
    storage_abrdf_index = obj.meta.abrdf_index_logical_to_storage(logical_abrdf_index);
    
    if storage_abrdf_index == -1
        chunk_index = -1;
        chunk_abrdf_index = -1;
        return;
    end

    chunk_index = floor((storage_abrdf_index - 1) / obj.meta.abrdfs_per_chunk) + 1;

    % determine offset inside current chunk
    chunk_abrdf_index = storage_abrdf_index - obj.meta.abrdfs_per_chunk * (chunk_index - 1);
end
