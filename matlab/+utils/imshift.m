% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-04-08
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
% Shift an image (i.e. a 2D or 3D array) along its first two dimensions by
% the given offsets. The data is wrapped around at the corners.
%
% Usage:
% image = imshift(image, offset_x, offset_y)
% 
% offset_x, offset_y are arbitrary pixel offsets that are cast to the
% nearest integers and wrap around at the image borders.
function img = imshift(img, offset_x, offset_y)
    
    h = size(img, 1);
    w = size(img, 2);
    
    xinds = int32(1 : w);
    yinds = int32(1 : h);
    xinds = mod(xinds - int32(offset_x) - 1, w) + 1;
    yinds = mod(yinds - int32(offset_y) - 1, h) + 1;
    img = img(yinds, xinds, :);
end
