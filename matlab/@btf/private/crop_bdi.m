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
% Crop a bidirectional image file (BDI).
function obj = crop_bdi(obj, roi, strides, output_file_name)
    meta_old = obj.meta;
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    width = obj.meta.width;
    mask = obj.meta.mask;
    
    if strcmpi(obj.meta.file_name, output_file_name)
        error('Cannot write cropped BDI to the same file %s', output_file_name);
    end
    
    try
        chunk_scan_lines = obj.meta.scan_lines_per_chunk;

        % x- and y-indices in the input BDI
        xs_in = roi(1, 1) : strides(1) : roi(1, 2);
        ys_in = roi(2, 1) : strides(2) : roi(2, 2);

        % texture size of the ROI and the output BDI
        roi_width = numel(xs_in);
        roi_height = numel(ys_in);

        if obj.meta.num_rotations > 0
            error('cropping BDIs with rotations is not implemented yet');
        end
        
        % update some of the meta data as we need those for writing the
        % BDI header
        obj.meta.width = roi_width;
        obj.meta.height = roi_height;
        obj.meta.mask = true(roi_width, roi_height);

        fid_cropped = fopen(output_file_name, 'w');
        obj.write_bdi_header(fid_cropped);

        % chunk information for cropped BDI
        chunk_begin = ftell(fid_cropped);
        num_stored_abrdfs = roi_width * roi_height;
        abrdfs_per_chunk = roi_width * chunk_scan_lines;
        num_chunks = ceil(num_stored_abrdfs / abrdfs_per_chunk);
        chunk_offset = zeros(num_chunks, 1, 'int64');

        if ~all(all(mask(ys_in, xs_in)))
            % the slow & painful way, where we have to account for missing
            % ABRDFs
            error('Cropping BDIs with masked ABRDFs is not implemented yet. All ABRDFs inside the ROI must be present.');
        else
            % chunk indices in the input BDI
            chunks_in = floor((ys_in - 1) / chunk_scan_lines) + 1;
            % scan line indices inside the chunks
            chunk_in_scan_lines =  mod(ys_in - 1, 4) + 1;
            % we allocate a buffer that can store two (cropped) chunks at a time
            chunk_out_scan_lines = mod((1 : roi_height) - 1, 2 * chunk_scan_lines) + 1;
            chunks_temp = zeros(nC, nL, nV, roi_width, 2 * chunk_scan_lines, 'uint16');
            
            % last output scan line that has been buffered
            max_scan_line_out = 0;
            % index of output chunk that has been written
            chunk_out_counter = 1;
            
            % iterate over input chunks
            chunks_in_unique = unique(chunks_in);
            obj.data.fid = fopen(obj.meta.file_name, 'r');
            for ci = 1 : numel(chunks_in_unique)
                obj.progress(ci / numel(chunks_in_unique), 'cropping BDI');
                chunk_in_index = chunks_in_unique(ci);
                
                % read input chunk
                chunk = reshape(obj.get_bdi_chunk(chunk_in_index, [], [], false), nC, nL, nV, width, []);

                % determine position where to put this in the output buffer
                out_scan_lines = chunk_out_scan_lines(chunks_in == chunk_in_index);

                % grab ROI and store it in output buffer
                chunks_temp(:, :, :, :, out_scan_lines) = ...
                    chunk(:, :, :, xs_in, chunk_in_scan_lines(chunks_in == chunk_in_index));

                % store last buffered scan line (special case at the end of ROI)
                max_scan_line_out = mod(max(max_scan_line_out, max(out_scan_lines)) - 1, chunk_scan_lines * 2) + 1;

                % if buffer is full, write it to file as two new chunks
                if max(out_scan_lines) == 2 * chunk_scan_lines || ci == numel(chunks_in_unique)
                    % first chunk
                    chunk_offset(chunk_out_counter) = ftell(fid_cropped);
                    chunk_out_counter = chunk_out_counter + 1;
                    fwrite(fid_cropped, numel(chunks_temp(:, :, :, :, 1 : chunk_scan_lines))  * utils.sizeof('uint16'), 'uint32');
                    fwrite(fid_cropped, typecast(reshape(chunks_temp(:, :, :, :, 1 : chunk_scan_lines), [], 1), 'uint8'), 'uint8');

                    % second chunk (can be less than chunk_scan_lines at the
                    % bottom of ROI)!
                    chunk_offset(chunk_out_counter) = ftell(fid_cropped);
                    chunk_out_counter = chunk_out_counter + 1;
                    fwrite(fid_cropped, numel(chunks_temp(:, :, :, :, chunk_scan_lines + 1 : max_scan_line_out))  * utils.sizeof('uint16'), 'uint32');
                    fwrite(fid_cropped, typecast(reshape(chunks_temp(:, :, :, :, chunk_scan_lines + 1 : max_scan_line_out), [], 1), 'uint8'), 'uint8');

                    max_scan_line_out = 0;
                end
            end

            obj.progress();
        end
        fclose(obj.data.fid);
        fclose(fid_cropped);
        
        % only now can we update the chunk meta information as we no longer
        % need the old one for reading the input chunks
        obj.meta.chunk_begin = chunk_begin;
        obj.meta.abrdfs_per_chunk = abrdfs_per_chunk;
        obj.meta.chunk_size = obj.meta.num_dirs * obj.meta.num_channels * abrdfs_per_chunk;
        obj.meta.num_stored_abrdfs = num_stored_abrdfs;
        obj.meta.abrdf_index_logical_to_storage = 1 : num_stored_abrdfs;
        obj.meta.num_chunks = num_chunks;
        obj.meta.chunk_offset = chunk_offset;
        obj.meta.file_name = output_file_name;
        
        obj.clear_buffer();
    catch err
        % ensure object is in a valid state in case of errors
        obj.meta = meta_old;
        error(err);
    end
end
