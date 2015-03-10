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
% Read uncompressed bidirectional image files (BDI) in Bonn University's binary
% format. A BDI file essentially stores the full BTF tensor in an ABRDF-wise
% manner, along with meta data like the bidirectional sampling. Multiple ABRDFs
% are put together into chunks. A chunk stores the ARBDFs corresponding to
% several scan lines in texture space. To access ABRDFs, we only need to
% determine the chunk index and the offset inside the chunk. Extracting textures
% is extremely inefficient, as we need to scan through all chunks, essentially
% reading the whole file. To speed things up, parts or all of the chunks can be
% buffered in memory using the method buffer_bdi(). This method attempts to
% store as many chunks contiguously in memory as possible, allowing for much
% faster access.
function [data, meta] = read_ubo_bdi(fid, signature, header_flag)
    fseek(fid, numel(signature), 'bof');
    
    meta = read_meta_data(fid);
    meta = read_bidir_sampling(fid, meta);
    meta.num_dirs = size(meta.L,1) * size(meta.V,1);
    meta.abrdf_size = meta.num_dirs * meta.num_channels; % number of elements in one ABRDF

    % read spatial extents
    meta.width = fread(fid, 1, 'uint32');
    meta.height = fread(fid, 1, 'uint32');

    % chunk info
    meta.scan_lines_per_chunk = fread(fid, 1, 'uint32');

    % there's a flag for the compression type
    meta.compression = 1; % 0: BDI_FULL_UNCOMPRESSED, 1: BDI_FULL_BZIP, 2: BDI_SPARSE_UNCOMPRESSED
    if header_flag
        meta.compression = fread(fid, 1, 'uchar');
    end
    
    if meta.compression ~= 2
        error('only BDIs of type BDI_SPARSE_UNCOMPRESSED are supported at the moment');
    end

    % this mask indicates which ABRDFs are missing
    meta.mask = fread_matrix(fid, 'uint8=>logical', meta.width, meta.height);

    % read rotations (if any)
    if strcmp(signature(end), 'R') || strcmp(signature(end - 2), 'R')
        meta.num_rotations = double(fread(fid, 1, 'uint32'));
        if (meta.num_rotations > 0)
            meta.rotations = fread_matrix(fid, 'float', meta.num_rotations, 3, 3);
        end
    else
        meta.num_rotations = 0;
    end

    % chunk info
    meta.chunk_begin = ftell(fid);
    meta.abrdfs_per_chunk = meta.width * meta.scan_lines_per_chunk;
    
    % total number of elements in each chunk
    meta.chunk_size = meta.num_dirs * meta.num_channels * meta.abrdfs_per_chunk;

    % BDI bzipped?
    meta.is_buffered = meta.compression == 1;

    % determine number of stored ABRDFs and their indices in case of a sparse
    % BDI
    if meta.compression == 2 % sparse uncompressed
        meta.num_stored_abrdfs = 0;
        meta.abrdf_index_logical_to_storage = zeros(numel(meta.mask), 1);
        for ai = 1 : numel(meta.mask)
            if meta.mask(ai) ~= 0
                meta.abrdf_index_logical_to_storage(ai) = meta.num_stored_abrdfs + 1;
                meta.num_stored_abrdfs = meta.num_stored_abrdfs + 1;
            else
                meta.abrdf_index_logical_to_storage(ai) = -1;
            end
        end
    else
        meta.num_stored_abrdfs = meta.width * meta.height;
        meta.abrdf_index_logical_to_storage = 1 : meta.num_stored_abrdfs;
    end

    % determine chunk offsets in file
    meta.num_chunks = ceil(meta.num_stored_abrdfs / meta.abrdfs_per_chunk);
    meta.chunk_offset = zeros(meta.num_chunks, 1, 'int64');

    % determine offsets of chunks
    max_compressed_chunk_size = 0;
    fseek(fid, meta.chunk_begin, 'bof');
    bytes_read_total = int64(meta.chunk_begin);
    for c = 1 : meta.num_chunks
        meta.chunk_offset(c) = bytes_read_total;
        compressed_chunk_size = fread(fid, 1, 'uint32');
        file_off = int64(compressed_chunk_size);
        fseek(fid, file_off, 'cof');
        if compressed_chunk_size > max_compressed_chunk_size
            max_compressed_chunk_size = compressed_chunk_size;
        end
        bytes_read_total = bytes_read_total + file_off + num_bytes('uint32');
    end

    % allocate storage for compressed chunks
    if meta.compression == 1 % full bzip compressed BDI
        meta.compressed_chunk = zeros(max_compressed_chunk_size, 1, 'char');
    end
    
    % this mask indicates which of the chunks have been buffered by buffer_bdi
    data.chunks_buffered = false(1, meta.num_chunks);
    data.num_chunks_in_buffer = 0;
    % disabling this leads to a hybrid mode where only those chunks that are not
    % buffered are read from file (this is still quite slow and therefore by
    % default only buffered data is displayed, whereas the rest is set black)
    data.only_use_buffered = true;
end
