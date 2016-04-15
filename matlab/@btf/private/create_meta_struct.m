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
function obj = create_meta_struct(obj, varargin)
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
    obj.meta.measurement_setup = '';
    obj.meta.image_sensor = '';
    obj.meta.light_source = '';
    obj.meta.ppmm = '';
    obj.meta.rgb_scale_factor = [1, 1, 1];
    obj.meta.xml = '';
    
    % set BTF / BDI dimensions
    obj.meta.nV = p.Results.num_views;
    obj.meta.nL = p.Results.num_lights;
    obj.meta.width = p.Results.width;
    obj.meta.height = p.Results.height;
    obj.meta.num_channels = 3;
    
    % channel and other information
    obj.meta.channel_names = 'RGB';
    obj.meta.num_rotations = 0;
    obj.meta.dynamic_range_reduction_method = p.Results.dynamic_range_reduction_method;
    
    % sampling information
    obj.meta.L = p.Results.view_dirs;
    obj.meta.V = p.Results.view_dirs;
    obj.meta.cosine_flag = p.Results.cosine_flag;
    
    % create default dome sampling if the input dimensions agree
    if isempty(obj.meta.L) || isempty(obj.meta.V)
        if obj.meta.nL == 151 && obj.meta.nV == 151
            [obj.meta.L, obj.meta.V] = utils.create_ubo_dome1_sampling();
        else
            error('please provide arrays for the angular sampling!');
        end
    end
    
    % handle all remaining meta data fields
    if ~isempty(fieldnames(p.Unmatched))
        for f = fieldnames(p.Unmatched)'
            obj.meta.(f{1}) = p.Unmatched.(f{1});
        end
    end
end
