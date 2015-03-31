% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-03-16
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
% Draw a rectangle into an image axes (or any other 2D plot).
%
% Usage:
% handles = imroi2(axes_handle, x_min, y_min, x_max, y_max, linestyle)
% 
% x_min, y_min, x_max, y_max: rectangle bounds in image coordinates
% linestyle: same syntax as understood by the plot command (e.g. 'r--' for a
% red, dotted line)
function handles = imroi2(axes_handle, x_min, y_min, x_max, y_max, linestyle)
    if ~exist('linestyle', 'var')
        linestyle = 'g-';
    end
    
    handles{1} = plot(axes_handle, ...
        [x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min], linestyle);
end
