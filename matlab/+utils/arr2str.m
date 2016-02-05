% *************************************************************************
% * Copyright 2016 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2016-02-05
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
% This function produces a nice looking string from a 1D or 2D array.
% Optionally, the output can be forced into one line. The created string can be
% copied and pasted to be used as matlab code again.
function str = arr2str(arr, one_line)
    if ~exist('one_line', 'var')
        one_line = true;
    end
    
    tmp = arrayfun(@(x) sprintf('%s, ', num2str(x)), arr, 'UniformOutput', false);
    tmp = mat2cell(tmp, ones(size(tmp, 1), 1), size(tmp, 2));
    tmp = cellfun(@(x) [x{:}], tmp, 'UniformOutput', false);
    
    if one_line
        tmp = cellfun(@(x) [x(1 : end - 2), '; '], tmp, 'UniformOutput', false);
    else
        tmp = cellfun(@(x) sprintf('%s\n', [x(1 : end - 2), ';']), tmp, 'UniformOutput', false);
    end
    
    tmp{end} = tmp{end}(1 : end - 2);
    str = ['[', tmp{:}, ']'];
end
