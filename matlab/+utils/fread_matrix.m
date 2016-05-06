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
% This method reads the first M rows and the first N columns of a matrix
% that is stored in row-major order in a file.
%
% mat = fread_matrix(fid, scalar_type, rows, columns, columns_file), where
% 
% - fid is the file id provided by fopen(), the current reading position
%   needs to be where the matrix starts
% - scalar_type is a string describing the data class as understood by
%   fread, this can also explicitly specify data type conversion, e.g.
%   uint8=>logical
% - rows and columns are the number of rows and columns that should be read
% - columns_file is the actual number of columns that are stored in the
%   file
function mat = fread_matrix(fid, scalar_type, rows, columns, columns_file)
    if ~exist('columns_file', 'var')
        columns_file = columns;
    end
    
    % handle cases like 'uint8=>logical', i.e. explicit data type conversion
    split_index = strfind(scalar_type, '=>');
    if ~isempty(split_index)
        scalar_type_mem = scalar_type(split_index + 2 : end);
        scalar_type_file = scalar_type(1 : split_index - 1);
    else
        % if only 'uint8' is specified to fread, an implicit conversion to
        % double occurs; to avoid this, we prepend '*'
        if scalar_type(1) ~= '*'
            scalar_type = ['*', scalar_type];
        end
        
        % conversion defaults to scalar type in file
        scalar_type_file = scalar_type(2 : end);
        scalar_type_mem = scalar_type_file;
    end

    num_skipped_cols = columns_file - columns;
    if num_skipped_cols > 0
        % this performs TERRIBLY poor but is left for better readability
%         % we're reading only parts of the matrix -> read row by row
%         mat = zeros(rows, columns, scalar_type(2:end));
%         for ri = 1 : rows
%             mat(ri, :) = fread(fid, columns, scalar_type);
%             fseek(fid, num_skipped_cols * size_of_scalar, 'cof');
%         end

        % get byte size of scalars in file
        size_of_scalar = utils.sizeof(scalar_type_file);
        
        mat = fread(fid, [rows, columns], ...
            sprintf('%d*%s=>%s', columns, scalar_type_file, scalar_type_mem), ...
            num_skipped_cols * size_of_scalar);
        mat = reshape(mat, columns, rows)';
    else
        % read the whole matrix at once
        % file stores row-major -> transpose after reading
        mat = reshape(fread(fid, rows * columns, scalar_type), columns, rows)';
    end
end
