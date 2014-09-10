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
% Write meta data information of a BTF file in one of Bonn University's
% binary formats.
function write_meta_data(fid, meta)
    pos_start = ftell(fid);
    
    struct_version = uint32(4);
    
    % we write three fixed length (80) strings
    measurement_setup = repmat(sprintf('\00'), 1, 80);
    image_sensor = repmat(sprintf('\00'), 1, 80);
    light_source = repmat(sprintf('\00'), 1, 80);
    measurement_setup(1 : numel(meta.measurement_setup)) = meta.measurement_setup;
    image_sensor(1 : numel(meta.image_sensor)) = meta.image_sensor;
    light_source(1 : numel(meta.light_source)) = meta.light_source;
        
    % predict the number of bytes of the whole meta data chunk
    size_bytes = sizeof('uint32') + sizeof(struct_version) + ...
        numel(meta.measurement_setup) + ...
        numel(meta.image_sensor) + ...
        numel(meta.light_source) + ...
        sizeof(meta.ppmm) + ...
        sizeof(meta.rgb_scale_factor);
    
    if struct_version > 1
        size_bytes = size_bytes + sizeof(meta.cosine_flag);
    end
        
    if struct_version > 2
        size_bytes = size_bytes + sizeof('uint32') + numel(meta.xml);
    end
    
    if struct_version > 3
        size_bytes = size_bytes + sizeof('uint32');
        for c = 1 : numel(meta.channel_names)
            size_bytes = size_bytes + sizeof('uint32') + numel(meta.channel_names{c});
        end
    end
    
    % start writing the actual meta data
    fwrite(fid, size_bytes, 'uint32');
    fwrite(fid, struct_version, 'uint32');
    fwrite(fid, measurement_setup, 'char');
    fwrite(fid, image_sensor, 'char');
    fwrite(fid, light_source, 'char');
    fwrite(fid, meta.ppmm, 'single');
    fwrite(fid, meta.rgb_scale_factor, 'single');
    
    % is the cosine of the light source's polar angle already multiplied
    % into the data?
    if struct_version > 1
        fwrite(fid, meta.cosine_flag, 'uint32');
    end
    
    % let's also store the big measurement XML string!
    if struct_version > 2
        str_size = numel(meta.xml);
        fwrite(fid, str_size, 'uint32');
        fwrite(fid, meta.xml, 'char');
    end
    
    % store color channel descriptors (this is important for multispectral
    % data!)
    if struct_version > 3
        fwrite(fid, meta.num_channels, 'uint32');
        for c = 1 : meta.num_channels
            str_size = numel(meta.channel_names{c});
            fwrite(fid, str_size, 'uint32');
            fwrite(fid, meta.channel_names{c}, 'char');
        end
    end
    
    % sanity check if we have actually written the expected number of bytes
    pos_end = ftell(fid);
    if pos_end - pos_start ~= size_bytes
        error('write_meta_data: size of written bytes not as expected (written: %d bytes, expected: %d bytes)', pos_end - pos_start, size_bytes);
    end
end