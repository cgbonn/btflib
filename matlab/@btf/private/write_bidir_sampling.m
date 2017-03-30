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
% Write bidirectional sampling to a file in one of Bonn University's binary
% formats.
function write_bidir_sampling(fid, meta)
    % so far only this is supported
    is_tensor_product_sampling = true;
    
    if is_tensor_product_sampling
        fwrite(fid, meta.nV, 'uint32');
        
        V_sph = utils.cart2sph2(meta.V(:, 1), meta.V(:, 2), meta.V(:, 3));
        L_sph = utils.cart2sph2(meta.L(:, 1), meta.L(:, 2), meta.L(:, 3))';
        for v = 1 : meta.nV
            fwrite(fid, V_sph(v, :), 'single');
            fwrite(fid, meta.nL, 'uint32');
            fwrite(fid, L_sph, 'single');
        end
    else
        error('write_bidir_sampling: not implemented for non-tensor-product samplings');
    end
end
