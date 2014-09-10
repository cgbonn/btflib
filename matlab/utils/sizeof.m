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
% Returns the number of bytes of an input array or of a scalar class, if
% the class name is provided as a string.
function n = sizeof(data)
    if ischar(data) && length(data) > 1
        if data(1) == '*'
            data = data(2 : end);
        end
        class_name = data;
        num_elements = 1;
    else
        class_name = class(data);
        num_elements = numel(data);
    end
    
    switch class_name
        case {'logical','uint8','int8'}
            n = 1;
        case {'char','uint16','int16'}
            n = 2;
        case {'uint32','int32'}
            n = 4;
        case {'uint64','int64'}
            n = 8;
        case 'single'
            n = 4;
        case 'double'
            n = 8;
        otherwise
            error('sizeof: unknown scalar type %s!', class_name)
    end
    n = n * num_elements;
end