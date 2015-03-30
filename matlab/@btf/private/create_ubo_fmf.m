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
% Given meta and data structs, this method constructs a FMF btf object. See the
% documentation of the btf constructor for further instructions.
function obj = create_ubo_fmf(obj, meta, data)
    obj.create_meta_struct(meta);
    
    w = obj.meta.width;
    h = obj.meta.height;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    nC = obj.meta.num_channels;
    
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('U', [], @(x) size(x, 1) == nC * w * h || size(x, 1) == nC * nL * nV);
    p.addParameter('SxV', [], @(x) size(x, 1) == w * h || size(x, 1) == nL * nV);
    p.addParameter('S', [], @(x) isnumeric(x));
    p.parse(data);
    
    obj.meta.num_components = size(p.Results.U, 2);
    
    % assign compontens
    if isempty(p.Results.S)
        obj.data.S = ones(obj.meta.num_components, 1);
    else
        obj.data.S = p.Results.S;
    end
    obj.data.U = p.Results.U;
    obj.data.SxV = p.Results.SxV;
    obj.data.class = class(obj.data.U);
end
