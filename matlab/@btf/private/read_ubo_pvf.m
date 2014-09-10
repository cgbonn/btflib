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
% * the terms of the GNU General Public License as published by the Free
% * Software Foundation, either version 3 of the License, or (at your
% * option) any later version.
% *
% * btflib is distributed in the hope that it will be useful, but WITHOUT
% * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% * for more details.
% *
% * You should have received a copy of the GNU General Public License along
% * with btflib.  If not, see <http://www.gnu.org/licenses/>.
% *
% *************************************************************************
%
% Read per view factorized (PVF) BTF file in Bonn University's binary
% format. In this format, for each view direction a PCA factorization is
% stored.
function [data, meta] = read_ubo_pvf(fid, signature, quality)
    fseek(fid, numel(signature), 'bof');
    magic_number = fread(fid, 1, 'int32');
    old_format = magic_number == 301074;

    if old_format
        % support for legacy PVF format
        num_components_in_file = fread(fid, 1, 'uint32');
        meta.width = double(fread(fid, 1, 'uint32'));
        meta.height = double(fread(fid, 1, 'uint32'));
        num_view_slots = fread(fid, 1, 'uint32');
        meta = read_bidir_sampling(fid, meta);
    else
        % currently used PVF format
        fseek(fid, -sizeof('int32'), 'cof');
        meta = read_common_header(fid, signature);
        num_components_in_file = fread(fid, 1, 'uint32');
        num_view_slots = fread(fid, 1, 'uint32');
    end

    data.EVs = cell(num_view_slots,1);
    data.Ms = cell(num_view_slots,1);
    data.Cs = cell(num_view_slots,1);
    data.Ws = cell(num_view_slots,1);
    for v = 1 : num_view_slots
        scalar_size = fread(fid, 1, 'char');
        data_type = size_to_type(scalar_size);

        num_components_file = fread(fid, 1, 'uint32');
        num_weights = fread(fid, 1, 'uint32');
        dimension = fread(fid, 1, 'uint32');

        % apply quality parameter by determining the actual number of
        % components to load
        num_components_to_load = max(1, round(min(1, quality) * num_components_file));
        
        % only read at most the number of components stored in file
        num_components = max(1, min(num_components_file, num_components_to_load));
        num_components_not_loaded = num_components_file - num_components;

        % read PCA eigenvalues, the subtracted mean, principal components and weights
        data.EVs{v} = fread_scalars(fid, num_components, data_type);
        fseek(fid, num_components_not_loaded * scalar_size, 'cof');
        data.Ms{v} = fread_scalars(fid, dimension, data_type);
        data.Cs{v} = fread_matrix(fid, data_type, dimension, num_components, num_components_file);
        data.Ws{v} = fread_matrix(fid, data_type, num_weights, num_components, num_components_file);
        
        % convert from half precision floats if necessary
        if scalar_size == 2
            data.EVs{v} = halfprecision(data.EVs{v}, 'single');
            data.Ms{v} = halfprecision(data.Ms{v}, 'single');
            data.Cs{v} = halfprecision(data.Cs{v}, 'single');
            data.Ws{v} = halfprecision(data.Ws{v}, 'single');
        end
    end
    if strcmpi(data_type, 'uint16')
        data.class = 'single';
    else
        data.class = data_type;
    end
end