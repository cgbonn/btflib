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
% Decode a texture from a per view factorized (PVF) BTF file in Bonn
% University's binary format. In this format, for each view direction a PCA
% factorization is stored.
function img = decode_pvf_texture(obj, l, v)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    h = obj.meta.height;
    w = obj.meta.width;
    Ms = obj.data.Ms;
    Cs = obj.data.Cs;
    Ws = obj.data.Ws;
    
    % determine transposition of BTF-matrix
    if size(Cs{1}, 1) == nC * nL
        % light fields stacked column-wise
        clInds = sub2ind([nC, nL], 1 : nC, repmat(l, 1, nC));
        img = permute(reshape(Ws{v} * Cs{v}(clInds, :)' + Ms{v}(clInds)', [w, h, nC]), [2, 1, 3]);
    elseif size(Cs{1}, 1) == nC * h * w
        % images stacked column-wise
        img = permute(reshape(Cs{v} * Ws{v}(l, :)' + Ms{v}, [nC, w, h]), [3,2,1]);
    else
        error('unknown format of BTF U-component');
    end
end