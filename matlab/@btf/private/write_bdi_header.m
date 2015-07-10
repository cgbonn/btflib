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
function write_bdi_header(obj, fid)
    % some parameters
    if isfield(obj.meta, 'num_scan_lines_per_chunk')
        num_scan_lines_per_chunk = obj.meta.num_scan_lines_per_chunk;
    else
        num_scan_lines_per_chunk = 4;
    end
    if isfield(obj.meta, 'compression')
        compression = obj.meta.compression;
    else
        compression = 2; % sparse uncompressed
    end
    if compression ~= 2
        error('only sparse uncompressed BDIs are supported at the moment.');
    end
    
    % write signature, meta data & bidir sampling
    fwrite(fid, '!BDIF06R2!', 'char');
    
    write_meta_data(fid, obj.meta);
    write_bidir_sampling(fid, obj.meta);
    
    % store spatial extents of textures (in pixels)
    fwrite(fid, obj.meta.width, 'uint32');
    fwrite(fid, obj.meta.height, 'uint32');
    
    fwrite(fid, num_scan_lines_per_chunk, 'uint32');
    if compression ~= 1
        fwrite(fid, compression, 'uint8');
    end
    
    % write sparsity mask (indicates which ABRDFs are missing)
    if ~isfield(obj.meta, 'mask')
        obj.meta.mask = true(obj.meta.width, obj.meta.height);
    end
    fwrite(fid, obj.meta.mask, 'uint8');
    
    % write rotations if available
    fwrite(fid, obj.meta.num_rotations, 'uint32');
    if obj.meta.num_rotations
        fwrite(fid, obj.meta.rotations', 'single');
    end
end
