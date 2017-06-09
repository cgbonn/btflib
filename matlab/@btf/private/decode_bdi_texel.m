% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2016-04-05
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
% Decode a single texel from a bidirectional image file (BDI).
function texel = decode_bdi_texel(obj, x, y, l, v)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    nlv = numel(l);
    assert(nlv == numel(v));
    nxy = numel(x);
    assert(nxy == numel(y));
    
    assert(nxy == nlv);
    
    abrdf = obj.decode_bdi_abrdf(x, y);
    texel = zeros(1, nC, nlv, 'single');
    for ci = 1 : nC
        inds = sub2ind([nL, nV, nxy, nC], l(:), v(:), ...
            (1 : nxy)', ci * ones(nlv, 1));
        texel(1, ci, :) = abrdf(inds);
    end
    texel = permute(texel, [3, 1, 2]);
end
