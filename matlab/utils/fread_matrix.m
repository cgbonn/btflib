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
% This method reads the first M rows and the first N columns of a matrix
% that is stored in row-major order in a file.
%
% mat = fread_matrix(fid, scalar_type, rows, columns, columns_file), where
% 
% - fid is the file id provided by fopen(), the current reading position
%   needs to be where the matrix starts
% - scalar_type is a string describing the data class as understood by
%   fread
% - rows and columns are the number of rows and columns that should be read
% - columns_file is the actual number of columns that are stored in the
%   file
function mat = fread_matrix(fid, scalar_type, rows, columns, columns_file)
    if ~exist('columns_file', 'var')
        columns_file = columns;
    end
    
    if scalar_type(1) ~= '*'
        scalar_type = ['*', scalar_type];
    end

    size_of_scalar = sizeof(scalar_type);

    num_skipped_cols = columns_file - columns;
    if num_skipped_cols > 0
        % this performs TERRIBLY poor but is left for better readability
%         % we're reading only parts of the matrix -> read row by row
%         mat = zeros(rows, columns, scalar_type(2:end));
%         for ri = 1 : rows
%             mat(ri, :) = fread(fid, columns, scalar_type);
%             fseek(fid, num_skipped_cols * size_of_scalar, 'cof');
%         end
        
        mat = fread(fid, [rows, columns], ...
            sprintf('%d*%s=>%s', columns, scalar_type(2:end), scalar_type(2:end)), ...
            num_skipped_cols * size_of_scalar);
        mat = reshape(mat, columns, rows)';
    else
        % read the whole matrix at once
        % file stores row-major -> transpose after reading
        mat = reshape(fread(fid, rows * columns, scalar_type), columns, rows)';
    end
end