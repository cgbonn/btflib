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
% Decode an apparent BRDF decorrelated full matrix factorized (DFMF) BTF
% file in Bonn University's binary format. In this format, the BTF's color
% values are  decorrelated in one of many alternative formats (e.g. the YUV
% color space). Each of the resulting decorrelated color channels is then
% arranged as a huge matrix which is factorized with an SVD. With this
% method, better compression ratios can be obtained by storing fewer
% components for some of the channels.
function abrdf = decode_dfmf_abrdf(obj, x, y)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    U = obj.data.U;
    SxV = obj.data.SxV;
    
    abrdf = zeros(nL * nV, nC, obj.data.class);
    % determine transposition of BTF-matrix
    if size(U{1},1) == nL * nV
        % abrdfs stacked column-wise
        xyInd = sub2ind([w, h], x, y);
        for c = 1 : nC
            abrdf(:, c) = U{c} * SxV{c}(xyInd, :)';
        end
    elseif size(U{1},1) == h * w
        % images stacked column-wise
        xyInd = sub2ind([w, h], x, y);
        warning('this is untested code!');
        for c = 1 : nC
            abrdf(:,c) = SxV{c} * U{c}(xyInd, :)';
        end
    else
        error('unknown format of BTF U-component');
    end
    
    abrdf = reshape(abrdf, [nL, nV, nC]);
end