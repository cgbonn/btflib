% *************************************************************************
% * Copyright 2015 University of Bonn
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
% Decode a single ABRDF from a bidirectional image file (BDI).
function abrdf = decode_bdi_abrdf(obj, x, y)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    n = numel(x);
    assert(numel(y) == n);
    x = x(:)';
    y = y(:)';
    
    abrdf = zeros(nC, nL, nV, n, 'single');
    for ii = 1 : n
        % to make this the least painful w.r.t. file access, x coordinates
        % should be unrolled first, then y. this is due to the fact that
        % chunks store several scan lines, which are along the x-dimension
        % of the textures. so use ndgrid(xs, ys) and not meshgrid(xs, ys).
        xi = x(ii);
        yi = y(ii);
        [chunk_index, chunk_abrdf_index] = obj.get_bdi_chunk_index(xi, yi);

        if chunk_index == -1
            continue;
        end

        if obj.data.chunks_buffered(chunk_index)
            abrdf_tmp = obj.data.chunks((chunk_abrdf_index - 1) * obj.meta.abrdf_size + 1: ...
                chunk_abrdf_index * obj.meta.abrdf_size, chunk_index);
        else
            obj.data.fid = fopen(obj.meta.file_name, 'r');
            abrdf_tmp = obj.get_bdi_chunk([], xi, yi);
            fclose(obj.data.fid);
        end

        % convert half precision floats to single precision floats & rearrange data
        abrdf(:, :, :, ii) = reshape(halfprecision(abrdf_tmp, 'single'), nC, nL, nV);
    end
    
    abrdf = permute(abrdf, [2, 3, 4, 1]);
end
