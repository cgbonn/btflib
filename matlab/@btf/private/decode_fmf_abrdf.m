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
% Decode an apparent BRDF from a full matrix factorized (FMF) BTF file in
% Bonn University's binary format. In this format, the BTF data is arranged
% as a huge matrix which is factorized with an SVD.
function abrdf = decode_fmf_abrdf(obj, x, y)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    U = obj.data.U;
    SxV = obj.data.SxV;
    n = numel(x);
    assert(numel(y) == n);
    
    % determine transposition of BTF-matrix
    if size(U,1) == nC * nL * nV
        % abrdfs stacked column-wise
        xyInd = sub2ind([w, h], x, y);
        abrdf = permute(reshape(U * SxV(xyInd, :)', [nC, nL, nV, n]), [2, 3, 4, 1]);
    elseif size(U,1) == nC * h * w
        % images stacked column-wise
        cxyInds = sub2ind([nC, w, h], 1 : nC, ...
            repmat(x, 1, nC), repmat(y, 1, nC));
        abrdf = permute(reshape(SxV * U(cxyInds, :)', [nL, nV, nC, n]), [1, 2, 4, 1]);
    else
        error('unknown format of BTF U-component');
    end
end
