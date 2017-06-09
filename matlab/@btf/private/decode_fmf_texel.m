% *************************************************************************
% * Copyright 2014 University of Bonn
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
% Decode a single texel from a full matrix factorized (FMF) BTF file in
% Bonn University's binary format. In this format, the BTF data is arranged
% as a huge matrix which is factorized with an SVD.
function texel = decode_fmf_texel(obj, x, y, l, v)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    U = obj.data.U;
    SxV = obj.data.SxV;
    nxy = numel(x);
    nlv = numel(l);
    assert(nxy == numel(y) && nlv == numel(v));
    
    texel = zeros(nlv, 1, nC, obj.data.class);
    % determine transposition of BTF-matrix
    if size(U, 1) == nC * nL * nV
        % abrdfs stacked column-wise
        xyInd = sub2ind([w, h], x, y);
        for ci = 1 : nC
            clvInds = sub2ind([nC, nL, nV], repmat(ci, 1, nlv), l, v);
            texel(:, 1, ci) = sum(U(clvInds, :) .* SxV(xyInd, :), 2);
        end
    elseif size(U,1) == nC * h * w
        % images stacked column-wise
        lvInds = sub2ind([nL, nV], l, v);
        for ci = 1 : nC
            cxyInds = sub2ind([nC, w, h], repmat(ci, 1, nlv), x, y);
            texel(:, 1, ci) = sum(SxV(lvInds, :) .* U(cxyInds, :), 2);
        end
    else
        error('unknown format of BTF U-component');
    end
end
