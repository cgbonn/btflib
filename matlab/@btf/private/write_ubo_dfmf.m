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
% Write decorrelated full matrix factorized (DFMF) BTF file in Bonn
% University's binary format. In this format, the BTF's color values are
% decorrelated in one of many alternative formats (e.g. the YUV color
% space). Each of the resulting decorrelated color channels is then
% arranged as a huge matrix which is factorized with an SVD. With this
% method, better compression ratios can be obtained by storing fewer
% components for some of the channels.
function write_ubo_dfmf(fid, data, meta)
    % write signature, meta data & bidir sampling
    write_common_header(fid, '!DFMF08FCR', meta);
    
    num_components = uint32(0);
    for c = 1 : meta.num_channels
        num_components = max(num_components, size(data.U{c}, 2));
    end
    
    % write DFMF maximum number of components & decorrelation method
    fwrite(fid, num_components, 'uint32');
    fwrite(fid, meta.color_model, 'int32');
    fwrite(fid, meta.color_mean, 'single');
    fwrite(fid, meta.color_transformation_matrix', 'single');
    
    % write decorrelated channels
    for c = 1 : meta.num_channels
        scalar_size = 2; % write as half precision floats
        data_type = size_to_type(scalar_size);

        fwrite(fid, scalar_size, 'uint8');
        fwrite(fid, size(data.U{c}, 2), 'uint32');
        fwrite(fid, size(data.U{c}, 1), 'uint32');
        fwrite(fid, size(data.SxV{c}, 1), 'uint32');

        % write the actual FMF matrices
        fwrite(fid, halfprecision(data.S{c}), data_type);
        fwrite(fid, halfprecision(data.U{c}'), data_type);
        fwrite(fid, halfprecision(data.SxV{c}'), data_type);
    end
end