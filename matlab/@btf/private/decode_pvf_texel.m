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
    nxy = numel(x);
    nlv = numel(l);
    assert(nxy == numel(y) && nlv == numel(v));
    
    % determine transposition of BTF-matrix
    if size(Cs{1}, 1) == nC * nL
        % light fields stacked column-wise
        texel = zeros(nxy, nc, nlv, obj.data.class);
        xyInd = sub2ind([h, w], y, x);
        for ii = 1 : nlv
            clInds = sub2ind([nC, nL], repmat((1 : nC)', 1, n), repmat(l(ii), 1, nC));
            texel(:, :, ii) = reshape(Ws{v(ii)}(xyInd, :) * Cs{v(ii)}(clInds, :) + ...
                Ms{v(ii)}(clInds)', nC, nxy);
        end
        texel = permute(texel, [3, 1, 2]);
    elseif size(Cs{1}, 1) == nC * h * w
        % images stacked column-wise
        texel = zeros(nC, nxy, nlv, obj.data.class);
        cxyInds = sub2ind([nC, h, w], repmat((1 : nC)', 1, nxy), ...
            repmat(y(:)', nC, 1), repmat(x(:)', nC, 1));
        for ii = 1 : nlv
            texel(:, :, ii) = reshape(Cs{v(ii)}(cxyInds, :) * ...
                Ws{v(ii)}(l(ii), :)', nC, nxy, 1) + Ms{v(ii)}(cxyInds);
        end
        texel = permute(texel, [3, 2, 1]);
    else
        error('unknown format of BTF U-component');
    end
end
