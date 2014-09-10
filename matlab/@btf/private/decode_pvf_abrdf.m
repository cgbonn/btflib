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
% Decode an apparent BRDF from a per view factorized (PVF) BTF file in Bonn
% University's binary format. In this format, for each view direction a PCA
% factorization is stored.
function abrdf = decode_pvf_abrdf(obj, x, y)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    Ms = obj.data.Ms; % means
    Cs = obj.data.Cs; % coefficients
    Ws = obj.data.Ws; % weights
    
    abrdf = zeros(nC * nL, nV, obj.data.class);
    % determine transposition of BTF-matrix
    if size(Cs{1}, 1) == nC * nL
        % light fields stacked column-wise
        xyInd = sub2ind([w, h], x, y);
        for v = 1 : nV
            abrdf(:, v) = Cs{v} * Ws{v}(xyInd, :)' + Ms{v};
        end
        abrdf = permute(reshape(abrdf, [nC, nL, nV]), [2, 3, 1]);
    elseif size(Cs{1}, 1) == nC * h * w
        % images stacked column-wise
        cxyInds = sub2ind([nC, w, h], 1 : nC, ...
            repmat(x, 1, nC), repmat(y, 1, nC));
        for v = 1 : nV
            abrdf(:, v) = reshape(Cs{v}(cxyInds, :) * Ws{v}' + ...
                repmat(Ms{v}(cxyInds), [1, nL]), [nC * nL, 1]);
        end
        abrdf = permute(reshape(abrdf, [nC, nL, nV]), [2, 3, 1]);
    else
        error('unknown format of BTF U-component');
    end
end