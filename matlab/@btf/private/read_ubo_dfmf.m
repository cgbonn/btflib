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
% Read decorrelated full matrix factorized (DFMF) BTF file in Bonn
% University's binary format. In this format, the BTF's color values are
% decorrelated in one of many alternative formats (e.g. the YUV color
% space). Each of the resulting decorrelated color channels is then
% arranged as a huge matrix which is factorized with an SVD. With this
% method, better compression ratios can be obtained by storing fewer
% components for some of the channels.
function [data, meta] = read_ubo_dfmf(fid, signature, quality)
    % decorrelated full matrix factorization
    fseek(fid, numel(signature), 'bof');
    meta = read_common_header(fid, signature);

    % read DFMF specific data
    meta.num_components = double(fread(fid, 1, 'uint32'));
    meta.color_model = fread(fid, 1, 'int32');
    meta.color_mean = fread(fid, 3, 'single');
    meta.color_transformation_matrix = fread_matrix(fid, 'single', 3, 3);

    % read compontens
    data.S = cell(meta.num_channels, 1);
    data.U = cell(meta.num_channels, 1);
    data.SxV = cell(meta.num_channels, 1);
    for c = 1 : meta.num_channels
        scalar_size = fread(fid, 1, 'uint8');
        data_type = size_to_type(scalar_size);

        % apply quality parameter by determining the actual number of
        % components to read
        num_components_file = fread(fid, 1, 'uint32');
        num_components = max(1, round(min(1, quality) * num_components_file));
        num_components_not_loaded = num_components_file - num_components;
        
        num_rows = fread(fid, 1, 'uint32');
        num_cols = fread(fid, 1, 'uint32');
        
        % singular values
        data.S{c} = fread_scalars(fid, num_components, data_type);
        fseek(fid, num_components_not_loaded * scalar_size, 'cof');
        
        % left singular vectors
        data.U{c} = fread_matrix(fid, data_type, num_rows, num_components, num_components_file);
        
        % right singular vectors (multiplied with singular values)
        data.SxV{c} = fread_matrix(fid, data_type, num_cols, num_components, num_components_file);
        
        % convert from half precision floats if necessary
        if (scalar_size == 2)
           data.S{c} = halfprecision(data.S{c}, 'single');
           data.U{c} = halfprecision(data.U{c}, 'single');
           data.SxV{c} = halfprecision(data.SxV{c}, 'single');
        end
    end
    
    if strcmpi(data_type, 'uint16')
        data.class = 'single';
    else
        data.class = data_type;
    end
end