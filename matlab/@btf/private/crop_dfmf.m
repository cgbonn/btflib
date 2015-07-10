% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-07-10
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
% Crop a decorrelated full matrix factorized (DFMF) BTF file in Bonn
% University's binary format. In this format, the BTF's color values are
% decorrelated in one of many alternative formats (e.g. the YUV color space).
% Each of the resulting decorrelated color channels is then arranged as a huge
% matrix which is factorized with an SVD. With this method, better compression
% ratios can be obtained by storing fewer components for some of the channels.
function obj = crop_dfmf(obj, roi, strides)
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    
    inds_x = roi(1, 1) : strides(1) : roi(1, 2);
    inds_y = roi(2, 1) : strides(2) : roi(2, 2);
    
    [inds_x, inds_y] = meshgrid(inds_x, inds_y);
    linds_xy = sub2ind([obj.meta.width, obj.meta.height], inds_x', inds_y');
    
    % determine transposition of BTF-matrix
    if size(obj.data.U{1}, 1) == nL * nV
        % images stacked row-wise
        for c = 1 : obj.meta.num_channels
            obj.data.SxV{c} = obj.data.SxV{c}(linds_xy, :);
        end
    elseif size(obj.data.U, 1) == h * w
        % images stacked column-wise
        for c = 1 : obj.meta.num_channels
            obj.data.U{c} = obj.data.U{c}(linds_xy, :);
        end
    else
        error('unknown format of BTF U-component');
    end
end
