% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-08-29
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
% Decode all eigen-textures from a per view factorized (PVF) BTF file in Bonn
% University's binary format. In this format, for each view direction a PCA
% factorization is stored. This method simply returns either the entire PCA
% coefficient or weight matrix of the PCA decomposition, depending on how
% the per-view matrices are arranged.
function textures = get_pvf_eigen_textures(obj)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    Ms = obj.data.Ms;
    Cs = obj.data.Cs;
    Ws = obj.data.Ws;
    
    % determine transposition of BTF-matrix
    textures = cell(nV, 1);
    if size(Cs{1}, 1) == nC * nL
        % light fields stacked column-wise
        for v = 1 : nV
            textures{v} = permute(reshape([ones(size(Ws{v}, 1), 1), Ws{v}], ...
                w, h, nC, []), [2, 1, 3, 4]);
        end
    elseif size(Cs{1}, 1) == nC * h * w
        % images stacked column-wise
        for v = 1 : nV
            textures{v} = permute(reshape([Ms{v}, Cs{v}], nC, w, h, []), [3, 2, 1, 4]);
        end
    else
        error('unknown format of BTF U-component');
    end
end