% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-03-10
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
% Read meta data information of a BTF file in one of Bonn University's
% binary formats.
function meta = read_meta_data(fid)
    pos_start = ftell(fid);

    % read meta data
    struct_size = fread(fid, 1, 'uint32'); % header length
    struct_version = fread(fid, 1, 'uint32');
    meta.measurement_setup = fread(fid, 80, '*char');
    meta.image_sensor = fread(fid, 80, '*char');
    meta.light_source = fread(fid, 80, '*char');
    meta.ppmm = fread(fid, 1, '*single');
    meta.rgb_scale_factor = fread(fid, 3, '*single');

    % is the cosine of the light source's polar angle already multiplied
    % into the data?
    if (struct_version > 1)
        meta.cosine_flag = fread(fid, 1, '*uint32');
    else
        meta.cosine_flag = false;
    end

    % all current formats also store the big measurement XML string
    if (struct_version > 2)
        xml_size = fread(fid, 1, '*uint32');
        if (xml_size > 0)
            meta.xml = fread(fid, xml_size, '*char');
        else
            meta.xml = '';
        end
        clear xml_size
    end

    % read color channel descriptors (this is important for multispectral
    % data!)
    meta.num_channels = 3;
    if (struct_version > 3)
        meta.num_channels = fread(fid, 1, 'uint32');
        meta.channel_names = repmat({''}, meta.num_channels, 1);
        for ci = 1 : meta.num_channels
            str_size = fread(fid, 1, '*uint32');
            meta.channel_names{ci} = fread(fid, str_size, '*char');
        end
        
        % deal with multispectral BTFs
        if meta.num_channels ~= 3
            try
                meta.wavelengths = str2double(meta.channel_names(:)');
                if any(isnan(meta.wavelengths))
                    error('non-numeric channel names?');
                end
            catch err
                error(['problem during conversion of channel names to double: %s', ...
                    '\nstored channel names are: ', repmat('''%s'' ', 1, numel(meta.channel_names))], ...
                    err.message, meta.channel_names{:});
            end
        end
    else
        meta.channel_names = 'RGB';
    end

    % sanity check if we have actually read the expected number of bytes
    pos_test = ftell(fid);
    if pos_start + struct_size ~= pos_test
        error('error reading meta data, was supposed to read %d bytes but actually read %d', struct_size, pos_test - pos_start);
    end
end
