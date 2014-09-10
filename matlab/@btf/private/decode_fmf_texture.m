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
% Decode a texture from a full matrix factorized (FMF) BTF file in Bonn
% University's binary format. In this format, the BTF data is arranged as a
% huge matrix which is factorized with an SVD.
function img = decode_fmf_texture(obj, l, v)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    U = obj.data.U;
    SxV = obj.data.SxV;
    
    % determine transposition of BTF-matrix
    if size(U, 1) == nC * nL * nV
        % abrdfs stacked column-wise
        clvInds = sub2ind([nC, nL, nV], 1 : nC, repmat(l, 1, nC), repmat(v, 1, nC));
        img = permute(reshape(SxV * U(clvInds, :)',[w, h, nC]), [2, 1, 3]);
    elseif size(U, 1) == nC * h * w
        % images stacked column-wise
        lvInd = sub2ind([nL, nV], l, v);
        img = permute(reshape(U * SxV(lvInd, :)', [nC, w, h]), [3, 2, 1]);
    else
        error('unknown format of BTF U-component');
    end
end