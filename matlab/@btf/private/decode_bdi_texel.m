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
    
    if ~obj.is_buffered() || ~all(obj.meta.mask(:))
        error('rendering from file is too slow, please buffer the entire BDI to memory');
    end
    
    if isa(obj.data.chunks, 'uint16')
        texel = zeros(nC, 1, nlv, 'uint16');
    else
        texel = zeros(nC, 1, nlv, class(obj.data.chunks));
    end
    for ci = 1 : nC
        inds = sub2ind([nC, nL, nV, obj.meta.width, ...
            obj.meta.scan_lines_per_chunk, obj.meta.num_chunks], ...
            repmat(ci, 1, nlv), l, v, x, ...
            mod(y - 1, obj.meta.scan_lines_per_chunk) + 1, ...
            floor((single(y) - 1) / obj.meta.scan_lines_per_chunk) + 1);
        texel(ci, 1, :) = obj.data.chunks(inds);
    end
    if isa(obj.data.chunks, 'uint16')
        texel = halfprecision(permute(texel, [3, 2, 1]), 'single');
    else
        texel = permute(texel, [3, 2, 1]);
    end
end
