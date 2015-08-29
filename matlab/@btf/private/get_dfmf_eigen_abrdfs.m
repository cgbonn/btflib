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
% Decode all eigen-ABRDFs from a decorrelated full matrix factorized
% (DFMF) BTF file in Bonn University's binary format. In this format, the
% BTF's color values are  decorrelated in one of many alternative formats
% (e.g. the YUV color space). Each of the resulting decorrelated color
% channels is then arranged as a huge matrix which is factorized with an
% SVD. With this method, better compression ratios can be obtained by
% storing fewer components for some of the channels. This method simply
% returns for each of the decorrelated channels either the entire U or V
% matrix of the SVD decomposition, depending on how the full matrix was
% arranged.
function abrdfs = get_dfmf_eigen_abrdfs(obj)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    U = obj.data.U;
    SxV = obj.data.SxV;
    
    % determine transposition of BTF-matrix
    if size(U{1}, 1) == nL * nV
        % imgs stacked column-wise
        abrdfs = U;
    elseif size(U, 1) == h * w
        % images stacked column-wise
        abrdfs = SxV;
    else
        error('unknown format of BTF U-component');
    end
    
    for c = 1 : nC
        abrdfs{c} = reshape(abrdfs{c}, nL, nV, []);
    end
end