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
% Given meta struct, this method sanity checks the struct values and assigns
% them to the meta struct of a btf object.
function meta = create_meta_struct(varargin)
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('width', -1, @(x) isscalar(x) && isnumeric(x) && x > 0);
    p.addParameter('height', -1, @(x) isscalar(x) && isnumeric(x) && x > 0);
    p.addParameter('num_lights', 151, @(x) isscalar(x) && isnumeric(x) && x > 0);
    p.addParameter('num_views', 151, @(x) isscalar(x) && isnumeric(x) && x > 0);
    p.addParameter('num_channels', 3, @(x) isscalar(x) && isnumeric(x) && x > 0);
    p.addParameter('light_dirs', [], @(x) isnumeric(x) && size(x, 2) == 3);
    p.addParameter('view_dirs', [], @(x) isnumeric(x) && size(x, 2) == 3);
    p.addParameter('cosine_flag', false, @(x) isscalar(x) && (isnumeric(x) || islogical(x)));
    p.addParameter('dynamic_range_reduction_method', 0, @(x) isscalar(x) && isnumeric(x));
    p.parse(varargin{:});
    
    % those are normally stored in the UBO BTF file formats
    meta.measurement_setup = '';
    meta.image_sensor = '';
    meta.light_source = '';
    meta.ppmm = 1;
    meta.rgb_scale_factor = [1, 1, 1];
    meta.xml = '';
    meta.file_name = '';
    
    % set BTF / BDI dimensions
    meta.nV = p.Results.num_views;
    meta.nL = p.Results.num_lights;
    meta.width = p.Results.width;
    meta.height = p.Results.height;
    meta.num_channels = p.Results.num_channels;
    
    % channel and other information
    meta.channel_names = {'R', 'G', 'B'};
    meta.num_rotations = 0;
    meta.dynamic_range_reduction_method = p.Results.dynamic_range_reduction_method;
    
    % sampling information
    meta.L = p.Results.light_dirs;
    meta.V = p.Results.view_dirs;
    
    if size(meta.L, 2) ~= 3
        warning('btf:light_transposed', 'light array should be nL x 3');
        meta.L = meta.L';
    end
    if size(meta.V, 2) ~= 3
        warning('btf:view_transposed', 'light array should be nV x 3');
        meta.V = meta.V';
    end
    assert(size(meta.L, 1) == meta.nL);
    assert(size(meta.V, 1) == meta.nV);
    
    meta.cosine_flag = p.Results.cosine_flag;
    
    % create default dome sampling if the input dimensions agree
    if isempty(meta.L) || isempty(meta.V)
        if meta.nL == 151 && meta.nV == 151
            [meta.L, meta.V] = utils.create_ubo_dome1_sampling();
        else
            error('please provide arrays for the angular sampling!');
        end
    end
    
    % handle all remaining meta data fields
    if ~isempty(fieldnames(p.Unmatched))
        for f = fieldnames(p.Unmatched)'
            meta.(f{1}) = p.Unmatched.(f{1});
        end
    end
    
    % ensure wavelengths are a row vector
    if isfield(p.Unmatched, 'wavelengths')
        wavelengths = p.Unmatched.wavelengths;
        if size(wavelengths, 1) ~= 1
            wavelengths = wavelengths';
        end
        meta.wavelengths = wavelengths;
    end
end
