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
% Replacement for Matlab's cart2sph with the following differences:
%   - instead of an elevation, this function computes a polar angle, i.e.
%     the angular distance in radians from the zenith (x = 0, y = 0, z = 1)
%   - per default, the radial coordinate is not computed as this function
%     is intended for normalized vectors; if it is required, an additional
%     boolean input parameter needs to be set to true
%   - input can either be a single N x 3 matrix which stores the cartesian
%     x-, y- and z-coordinates in its rows, or separate N x 1 arrays of
%     these three quantities
%   - it either returns an N x 2 array (default) with the polar angle in
%     its first and the azimuth angle in its second column; if requested,
%     it returns an N x 3 array which furthermore stores the radial
%     component in its third column
function sph = cart2sph2(varargin)
    % input parsing
    if numel(varargin) < 3
        % one input matrix of dimensions n x 3
        x = varargin{1}(:, 1);
        y = varargin{1}(:, 2);
        z = varargin{1}(:, 3);
    elseif numel(varargin) >= 3
        % three separate matrices
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
    end
    n = size(x, 1);
    if numel(varargin) > 3
        with_radius = varargin{4};
    else
        with_radius = false;
    end
    
    len_xy = hypot(x, y);
    
    % radius?
    if with_radius
        sph = zeros(n, 3, class(x));
        
        sph(:, 3) = hypot(len_xy, z);
    else
        sph = zeros(n, 2, class(x));
    end
    
    % inclination
    sph(:, 1) = atan2(len_xy, z);
    
    % azimuth
    sph(:, 2) = atan2(y, x);
end