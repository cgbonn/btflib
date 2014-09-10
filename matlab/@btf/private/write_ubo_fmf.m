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
% Write full matrix factorized (FMF) BTF file in Bonn University's binary
% format. In this format, the BTF data is arranged as a huge matrix which
% is factorized with an SVD.
function write_ubo_fmf(fid, data, meta)
    % write signature, meta data & bidir sampling
    write_common_header(fid, '!FMF12FCER', meta);
    
    header_version = char(2);
    
    fwrite(fid, header_version, 'uchar');
    fwrite(fid, meta.dynamic_range_reduction_method, 'uint32');
    
    % write height map if present
    if isfield(data, 'height_map') && ~isempty(data.height_map)
        fwrite(fid, [size(data.height_map, 2), size(data.height_map, 1)], 'uint32');
        fwrite(fid, halfprecision(data.height_map'), 'uint16');
    else
        fwrite(fid, [0,0], 'uint32');
    end
    
    % write FMF dimensions & data type
    num_components = size(data.U, 2);
    fwrite(fid, num_components, 'uint32');
    
    % write as half precision floats
    scalar_size = 2;
    data_type = size_to_type(scalar_size);
    
    fwrite(fid, scalar_size, 'uchar');
    fwrite(fid, num_components, 'uint32');
    fwrite(fid, size(data.U, 1), 'uint32');
    fwrite(fid, size(data.SxV, 1), 'uint32');
    
    % write the actual FMF matrices
    fwrite(fid, halfprecision(data.S), data_type);
    fwrite(fid, halfprecision(data.U'), data_type);
    fwrite(fid, halfprecision(data.SxV'), data_type);
end