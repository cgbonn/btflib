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
        texel = zeros(nlv, 1, nc, obj.data.class);
        xyInd = sub2ind([w, h], x, y);
        for ii = 1 : nlv
            for ci = 1 : nC
                clInd = sub2ind([nC, nL], ci, l(ii));
                texel(ii, 1, ci) = sum(Ws{v(ii)}(xyInd(ii), :) .* Cs{v(ii)}(clInd, :), 2) + ...
                    Ms{v(ii)}(clInd);
            end
        end
    elseif size(Cs{1}, 1) == nC * h * w
        % images stacked column-wise
        texel = zeros(nlv, 1, nC, obj.data.class);
        for ii = 1 : nlv
            for ci = 1 : nC
                cxyInd = sub2ind([nC, w, h], ci, x(ii), y(ci));
                texel(ii, 1, ci) = sum(Cs{v(ii)}(cxyInd, :) .* ...
                    Ws{v(ii)}(l(ii), :), 2) + Ms{v(ii)}(cxyInd);
            end
        end
        texel = permute(texel, [3, 2, 1]);
    else
        error('unknown format of BTF U-component');
    end
end
