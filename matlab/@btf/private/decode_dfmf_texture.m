% *************************************************************************
% * Copyright 2014 University of Bonn
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
% Decode a texture from a decorrelated full matrix factorized (DFMF) BTF
% file in Bonn University's binary format. In this format, the BTF's color
% values are  decorrelated in one of many alternative formats (e.g. the YUV
% color space). Each of the resulting decorrelated color channels is then
% arranged as a huge matrix which is factorized with an SVD. With this
% method, better compression ratios can be obtained by storing fewer
% components for some of the channels.
function img = decode_dfmf_texture(obj, l, v)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    U = obj.data.U;
    SxV = obj.data.SxV;
    
    % determine transposition of BTF-matrix
    img = zeros(h * w, nC, obj.data.class);
    if size(U{1}, 1) == nL * nV
        % abrdfs stacked column-wise
        lvInd = sub2ind([nL, nV], l, v);
        for c = 1 : nC
            img(:, c) = SxV{c} * U{c}(lvInd, :)';
        end
    elseif size(U{1}, 1) == h * w
        % images stacked column-wise
        lvInd = sub2ind([nL, nV], l, v);
        for c = 1 : nC
            img(:, c) = U{c} * SxV{c}(lvInd, :)';
        end
    else
        error('unknown format of BTF U-component');
    end
    
    img = permute(reshape(img, [w, h, nC]), [2, 1, 3]);
end
