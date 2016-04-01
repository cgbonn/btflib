% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-07-10
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
% Write uncompressed bidirectional image files (BDI) in Bonn University's binary
% format. A BDI file essentially stores the full BTF tensor in an ABRDF-wise
% manner, along with meta data like the bidirectional sampling. Multiple ABRDFs
% are put together into chunks. A chunk stores the ARBDFs corresponding to
% several scan lines in texture space.
% This function can in theory be used to write BDIs from non-BDI (i.e.
% compressed) BTFs, as long as those provide a decode_abrdf method.
function write_bdi_chunks(obj, fid)
    num_scan_lines_per_chunk = 4;
    abrdfs_per_chunk = obj.meta.width * num_scan_lines_per_chunk;
    abrdf_size = obj.meta.nL * obj.meta.nV * obj.meta.num_channels;
    
    logical_abrdf_indices = find(obj.meta.mask');
    num_chunks = ceil(numel(logical_abrdf_indices) / abrdfs_per_chunk);
    
    obj.data.fid = fopen(obj.meta.file_name, 'r');
    for chunk_index = 1 : num_chunks
        fprintf('extracting chunk %d / %d... ', chunk_index, num_chunks);
        
        if strcmp(obj.format_str, 'BDI')
            % get whole chunk at once
            if obj.data.chunks_buffered(chunk_index)
                chunk = obj.data.chunks(:, chunk_index);
            else
                chunk = obj.get_bdi_chunk(chunk_index, [], [], false);
            end
        else
            % fill chunk manually ABRDF-wise; this would potentially allow
            % writing BDIs from compressed BTFs
            chunk_start = (chunk_index - 1) * abrdfs_per_chunk + 1;
            chunk_end = min(numel(logical_abrdf_indices), chunk_index * abrdfs_per_chunk);
            abrdf_indices_chunk = logical_abrdf_indices(chunk_start : chunk_end)';
            chunk = zeros(abrdf_size, numel(abrdf_indices_chunk), 'uint16');
            for ai = 1 : numel(abrdf_indices_chunk)
                a = abrdf_indices_chunk(ai);
                [x, y] = ind2sub([obj.meta.width, obj.meta.height], a);
                abrdf = obj.decode_abrdf(x, y);
                abrdf = halfprecision(reshape(permute(abrdf, [3, 1, 2]), [], 1));
                chunk(:, ai) = abrdf;
            end
        end
        
        fprintf('writing to file...\n');
        fwrite(fid, numel(chunk)  * utils.sizeof('uint16'), 'uint32');
        % FIXME: why the typecasting?
        fwrite(fid, typecast(chunk(:), 'uint8'), 'uint8');
    end
    fclose(obj.data.fid);
end
