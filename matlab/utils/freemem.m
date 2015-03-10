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
% Platform independent function to query for the number of bytes of available
% memory.
function [mem_free,mem_total] = freemem()
    if isunix
        [~,w] = unix('free -b | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        mem_total = stats(1);
        mem_free = (stats(3)+stats(end));
    else
        [user,sys] = memory();
        mem_total = sys.PhysicalMemory.Total;
        mem_free = user.MaxPossibleArrayBytes;
    end
end
