% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-03-30
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
% Given meta and data structs, this method constructs a BDI btf object. See the
% documentation of the btf constructor for further instructions.
function obj = create_ubo_bdi(obj, meta, data)
    obj.create_meta_struct(meta);
    
    w = obj.meta.width;
    h = obj.meta.height;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    nC = obj.meta.num_channels;
    
    obj.meta.num_dirs = size(obj.meta.L, 1) * size(obj.meta.V, 1);
    obj.meta.abrdf_size = obj.meta.num_dirs * obj.meta.num_channels; % number of elements in one ABRDF

    % additional meta data
    obj.meta.compression = 2; % BDI_SPARSE_UNCOMPRESSED
    obj.meta.num_rotations = 0;
    
    % chunk info
    obj.meta.scan_lines_per_chunk = 4;
    obj.meta.abrdfs_per_chunk = obj.meta.width * obj.meta.scan_lines_per_chunk;
    obj.meta.chunk_begin = -1;
    obj.meta.chunk_size = obj.meta.num_dirs * obj.meta.num_channels * obj.meta.abrdfs_per_chunk;
    obj.meta.num_stored_abrdfs = obj.meta.width * obj.meta.height;
    obj.meta.abrdf_index_logical_to_storage = 1 : obj.meta.width * obj.meta.height;
    obj.meta.num_chunks = ceil(obj.meta.num_stored_abrdfs / obj.meta.abrdfs_per_chunk);
    obj.meta.chunk_offset = zeros(obj.meta.num_chunks, 1, 'int64');
    obj.meta.mask = true(obj.meta.width, obj.meta.height);
    
    % if a 5D BDI tensor is provided, reshape it to chunk format
    if isnumeric(data) && ndims(data) == 5
        if size(data, 1) ~= nC
            d = find(size(data) == nC, 1);
            data = permute(data, [d, setdiff(1 : 5, d)]);
        end
        data = reshape(data, nC * nL * nV, w * h);
    end
    
    % ensure data type is correct
    if isnumeric(data)
        switch class(data)
            case {'single', 'double'}
                data = halfprecision(data);
            case {'int16', 'uint16'}
                % half precision floats are stored as uint16 in matlab
            otherwise
                error(['cannot convert data to half precision float, ', ...
                    'data must be single or double precision floats, ', ...
                    'or already converted to half (i.e. int16 or uint16).']);
        end
    end
    
    p = inputParser;
    p.KeepUnmatched = true;
    p.addRequired('chunks', @(x) isnumeric(x) && ...
        (ismatrix(x) && all(size(x) == [nC * nL * nV, w * h]) || ...
        ismatrix(x) && all(size(x) == [obj.meta.chunk_size, obj.meta.num_chunks]) || ...
        ndims(x) == 5 && all(size(x) == [nC, nL, nV, w, h])));
    p.parse(data);
    
    obj.data.chunks = p.Results.chunks;
    clear p;
    
    if size(obj.data.chunks, 1) ~= obj.meta.chunk_size
        if numel(obj.data.chunks) ~= obj.meta.chunk_size * obj.meta.num_chunks
            obj.data.chunks = reshape(obj.data.chunks, [nC, nL, nV, w, h]);
            scan_lines_missing = mod(obj.meta.height, obj.meta.scan_lines_per_chunk);
            obj.data.chunks = cat(5, obj.data.chunks, ...
                zeros(nC, nL, nV, w, scan_lines_missing, 'like', obj.data.chunks));
        end
        obj.data.chunks = reshape(obj.data.chunks, obj.meta.chunk_size, obj.meta.num_chunks);
    end

    % this mask indicates which of the chunks have been buffered by buffer_bdi
    obj.data.chunks_buffered = true(1, obj.meta.num_chunks);
    obj.data.num_chunks_in_buffer = obj.meta.num_chunks;
    % diable loading textures from file
    obj.data.textures_from_file = false;
    obj.data.only_use_buffered = true;
end
