% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-03-31
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
% Draw a small crosshair onto an image axes (or any other 2D plot) and add text
% next  to the crosshair.
%
% Usage:
% handles = imannotate(axes_handle, x, y, linestyle, d1, d2, string, text_color)
% 
% linestyle:	same syntax as understood by the plot command (e.g. 'r--' for a
%               red, dotted line)
% d1, d2:       crosshair consists of four bars around the center coordinates,
%               d1 is the inner, d2 the outer distance to the selected
%               coordinates, i.e. the crosshair is in total 2 * d2 + 1 pixels in
%               diameter, with a 2 * d1 - 1 "hole" in the center
% text_color:   color as understood by the text function, i.e. a 3 element array
%               or a shorhand color string like 'r' for red text
function handles = imannotate(axes_handle, x, y, linestyle, d1, d2, string, text_color)
    if isempty(axes_handle)
        axes_handle = gca;
    end
    
    if ~exist('linestyle', 'var')
        linestyle = 'g-';
    end
    
    if ~exist('d1', 'var')
        d1 = 3;
    end
    
    if ~exist('d2', 'var')
        d2 = 5;
    end
    
    if ~exist('text_color', 'var')
        text_color = 'g';
    end
    
    handles{1} = plot(axes_handle, [x - d2, x - d1], [y, y], ...
        linestyle, 'LineWidth', 2);
    handles{2} = plot(axes_handle, [x + d1, x + d2], [y, y], ...
        linestyle, 'LineWidth', 2);
    handles{3} = plot(axes_handle, [x, x], [y - d2, y - d1], ...
        linestyle, 'LineWidth', 2);
    handles{4} = plot(axes_handle, [x, x], [y + d1, y + d2], ...
        linestyle, 'LineWidth', 2);
    
    if exist('string', 'var')
        handles{5} = text(x + d2 + 1, y + d2 + 1, string, ...
            'Parent', axes_handle, 'Color', text_color);
    end
end
