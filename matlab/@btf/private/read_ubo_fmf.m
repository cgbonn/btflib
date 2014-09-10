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
% Read full matrix factorized (FMF) BTF file in Bonn University's binary
% format. In this format, the BTF data is arranged as a huge matrix which
% is factorized with an SVD.
function [data, meta] = read_ubo_fmf(fid, signature, header_flag, quality)
    fseek(fid, numel(signature), 'bof');
    meta = read_common_header(fid, signature);

    header_version = 0;
    if header_flag
        header_version = fread(fid, 1, 'uchar');
    end

    % the dynamic range of the data might have been reduced to facilitate
    % the matrix factorization
    if header_version >= 1
        meta.dynamic_range_reduction_method = fread(fid, 1, 'uint32');
    end

    % the data could be resampled to the surface geometry, which we then
    % need to read as a height map
    if header_version >= 2
        image_size = fread(fid, 2, 'uint32');

        if prod(image_size) > 0
            data.height_map = fread_matrix(fid, 'uint16', image_size(1), image_size(2));
            data.height_map = halfprecision(data.height_map, 'single');
        end
    end

    if header_version >= 3
        error('unknown header version %d in DFMF', header_version);
    end

    % apply quality parameter -> determine actual number of components to load
    num_components_file = fread(fid, 1, 'uint32');
    num_components_to_load = max(1, round(min(1, quality) * num_components_file));

    % determine format of the parameters
    scalar_size = fread(fid, 1, 'char');
    data_type = size_to_type(scalar_size);
    if strcmpi(data_type, 'uint16')
        data.class = 'single';
    else
        data.class = data_type;
    end

    % double check if we're reading an allowed number of components
    num_components_file = fread(fid, 1, 'uint32');
    num_components = max(1, min(num_components_file, num_components_to_load));
    data.num_components = num_components;
    num_components_not_loaded = num_components_file - num_components;
    
    % dimensions of matrices in file
    num_rows = fread(fid, 1, 'uint32');
    num_cols = fread(fid, 1, 'uint32');

    % singular values
    data.S = fread_scalars(fid, num_components, data_type);
    fseek(fid, num_components_not_loaded * scalar_size, 'cof');

    % left singular vectors
    data.U = fread_matrix(fid, data_type, num_rows, num_components, num_components_file);

    % right singular vectors (multiplied by singular values)
    data.SxV = fread_matrix(fid, data_type, num_cols, num_components, num_components_file);
    
    % convert from half precision floats if necessary
    if scalar_size == 2
        data.S = halfprecision(data.S, 'single');
        data.U = halfprecision(data.U, 'single');
        data.SxV = halfprecision(data.SxV, 'single');
    end
end