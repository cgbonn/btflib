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
% Given meta and data structs, this method constructs a DFMF btf object. See the
% documentation of the btf constructor for further instructions.
function obj = create_ubo_dfmf(obj, meta, data)
    obj.create_meta_struct(meta);
    
    w = obj.meta.width;
    h = obj.meta.height;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    nC = obj.meta.num_channels;
    
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('color_model', 0, @(x) isnumeric(x) && isscalar(x));
    p.addParameter('color_mean', [0, 0, 0], @isnumeric);
    p.addParameter('color_transformation_matrix', eye(3), @isnumeric);
    p.parse(meta);
    
    % store information regarding color model
    obj.meta.color_model = p.Results.color_model;
    obj.meta.color_mean = p.Results.color_mean;
    obj.meta.color_transformation_matrix = p.Results.color_transformation_matrix;

    p = inputParser;
    p.addParameter('U', [], @(x) iscell(x) && numel(x) == nC && ...
        all(cellfun(@(y) size(y, 1) == w * h || size(y, 1) == nL * nV, x)));
    p.addParameter('SxV', [], @(x) iscell(x) && numel(x) == nC && ...
        all(cellfun(@(y) size(y, 1) == w * h || size(y, 1) == nL * nV, x)));
    p.addParameter('S', [], @(x) (isnumeric(x) && numel(x) == nC) || (iscell(x) ...
        && all(cellfun(@numel, x) == nC)));
    p.parse(data);
    
    % assign compontens
    obj.meta.num_components = size(p.Results.U{1}, 2);
    if isempty(p.Results.S)
        obj.data.S = repmat({ones(obj.meta.num_components, 1)}, nC, 1);
    else
        obj.data.S = p.Results.S;
    end
    obj.data.U = p.Results.U;
    obj.data.SxV = p.Results.SxV;
    obj.data.class = class(obj.data.U{1});
end
