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
% Decode a single ABRDF or buffer an entire chunk of ABRDFs from a bidirectional
% image file (BDI).
function abrdf = get_bdi_chunk(obj, chunk_index, x, y)
    if exist('x', 'var') && exist('y', 'var') && ~isempty(x) && ~isempty(y)
        % extract only a single ABRDF
        [chunk_index, chunk_abrdf_index] = obj.get_bdi_chunk_index(x, y);
        
        file_offset = obj.meta.chunk_offset(chunk_index);
        
        % abrdf is missing in sparse BDI
        if chunk_abrdf_index == -1
            abrdf = zeros(obj.meta.num_channels * obj.meta.nL * obj.meta.nV, 1, 'uint16');
            return;
        end
        
        if chunk_abrdf_index ~= 0
            file_offset = file_offset + 4;
        end
        
        % data is stored as half precision floats (which are stored as
        % uin16 in matlab)
        abrdf_byte_size = obj.meta.abrdf_size * sizeof('uint16');
        if chunk_abrdf_index == 0
            1;
        end
        file_offset = file_offset + (chunk_abrdf_index - 1) * abrdf_byte_size;
        
        read_size = obj.meta.abrdf_size;
    else
        % extract whole chunk
        file_offset = obj.meta.chunk_offset(chunk_index);
        file_offset = file_offset + sizeof('uint32');
        read_size = obj.meta.chunk_size;
    end
    
    fseek(obj.data.fid, file_offset, 'bof');
    
    % write the data directly into the buffer or return those chunks that don't
    % fit into memory
    if chunk_index < obj.data.num_chunks_in_buffer    
        obj.data.chunks(:, chunk_index) = fread(obj.data.fid, read_size, '*uint16', 0);
    elseif chunk_index == obj.data.num_chunks_in_buffer
        % the last chunk can be smaller than the rest, indexing into data.chunks
        % is slow, so we do it only here
        chunk = fread(obj.data.fid, read_size, '*uint16', 0);
        obj.data.chunks(1 : numel(chunk), chunk_index) = chunk;
    else
        % there was too little memory assigned to the chunk buffer, so we cannot
        % store chunks with indices larger than num_chunks_in_buffer
        abrdf = fread(obj.data.fid, read_size, '*uint16', 0);
    end
end
