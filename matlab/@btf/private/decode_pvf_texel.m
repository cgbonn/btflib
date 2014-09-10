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
% Decode a single texel from a per view factorized (PVF) BTF file in Bonn
% University's binary format. In this format, for each view direction a PCA
% factorization is stored.
function texel = decode_pvf_texel(obj, x, y, l, v)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    h = obj.meta.height;
    w = obj.meta.width;
    Ms = obj.data.Ms; % means
    Cs = obj.data.Cs; % coefficients
    Ws = obj.data.Ws; % weights
    
    % determine transposition of BTF-matrix
    if size(Cs{1}, 1) == nC * nL
        % light fields stacked column-wise
        clInds = sub2ind([nC, nL], 1 : nC, repmat(l, 1, nC));
        xyInd = sub2ind([h, w], y, x);
        texel = squeeze(Cs{v}(clInds, :) * Ws{v}(xyInd, :)' + Ms{v}(clInds));
    elseif size(Cs{1}, 1) == nC * h * w
        % images stacked column-wise
        cxyInds = sub2ind([nC, h, w], 1 : nC, ...
            repmat(y, 1, nC),repmat(x, 1, nC));
        texel = squeeze(Cs{v}(cxyInds, :) * Ws{v}(l, :)' + Ms{v}(cxyInds));
    else
        error('unknown format of BTF U-component');
    end
end