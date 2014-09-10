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
% Replacement for Matlab's sph2cart with the following differences:
%   - instead of an elevation, this function expects a polar angle, i.e.
%     the angular distance in radians from the zenith (x = 0, y = 0, z = 1)
%   - the radial coordinate is optional, if omitted, it is assumed to be 1
%   - input can either be a single N x 3 or N x 2 matrix which stores the
%     polar angle in its first, the azimuth angle in its second and
%     optionally the radius in its third column, or separate N x 1 arrays
%     of these three quantities
%   - it returns an N x 3 array of cartesian coordinates
%
% Usage:
%   xyz = sph2cart2(sph)
%   xyz = sph2cart2(polar, azimuth)
%   xyz = sph2cart2(polar, azimuth, radius)
% 
% where all input and out arrays are N x 3 or N x 1 respectively.
function xyz = sph2cart2(varargin)
    % input parsing
    if numel(varargin) == 1
        % one input matrix of dimensions n x 2 or n x 3
        inclination = varargin{1}(:, 1);
        n = size(inclination, 1);
        azimuth = varargin{1}(:, 2);
        
        if size(varargin{1}, 2) == 3
            radius = varargin{1}(:, 3);
        else
            radius = ones(n, 1, class(inclination));
        end
    elseif numel(varargin) >= 2
        % separate input matrices
        inclination = varargin{1};
        n = size(inclination, 1);
        azimuth = varargin{2};
        
        if numel(varargin) >= 3
            radius = varargin{3};
        else
            radius = ones(n, 1, class(inclination));
        end
    end
    
    xyz = zeros(n, 3, class(azimuth));
    
    rsinelev = radius .* sin(inclination);
    xyz(:, 1) = rsinelev .* cos(azimuth);
    xyz(:, 2) = rsinelev .* sin(azimuth);
    xyz(:, 3) = radius .* cos(inclination);
end