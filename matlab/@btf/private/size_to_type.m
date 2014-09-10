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
% Given the byte size of a scalar, returns the numeric class that is used
% in UBO btf files.
function data_type = size_to_type(scalar_size)
    if (scalar_size == 2)
        % this actually corresponds to half precision floats, which we need
        % to represent as uint16 in matlab
        data_type = 'uint16';
    elseif (scalar_size == 4)
        data_type = 'single';
    elseif (scalar_size == 8)
        data_type = 'double';
    else
        error('unknown datatype');
    end
end
