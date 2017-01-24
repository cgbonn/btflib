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
% Given meta and data structs, this method constructs a PVF btf object. See the
% documentation of the btf constructor for further instructions.
function obj = create_ubo_pvf(obj, meta, data)
    obj.meta = utils.create_meta_struct(meta);
    
    w = obj.meta.width;
    h = obj.meta.height;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    nC = obj.meta.num_channels;
    
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('Cs', [], @(x) iscell(x) && numel(x) == nV && ...
        all(cellfun(@(y) (size(y, 1) == nC * nL || size(y, 1) == nC * h * w), x)));
    p.addParameter('Ws', [], @(x) iscell(x) && numel(x) == nV && ...
        all(cellfun(@(y) size(y, 1) == w * h || size(y, 1) == nL, x)));
    p.addParameter('Ms', [], @(x) iscell(x) && ...
        all(cellfun(@(y) (size(y, 1) == nC * nL || size(y, 1) == nC * h * w), x)));
    p.addParameter('EVs', [], @(x) iscell(x));
    p.parse(data);

    % assign compontens
    obj.data.Ms = p.Results.Ms;
    obj.data.Cs = p.Results.Cs;
    obj.data.Ws = p.Results.Ws;
    if isempty(p.Results.EVs)
        obj.data.EVs = repmat({ones(size(obj.data.Cs{1}, 2))}, nV, 1);
    else
        obj.data.EVs = p.Results.EVs;
    end
    obj.data.class = class(obj.data.Cs{1});
end
