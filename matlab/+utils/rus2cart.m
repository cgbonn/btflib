% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-12-21
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
% This function computes the light and view directions in cartesian coordinates
% from half and difference vectors in the Rusinkiewicz BRDF parameterization.
% Input are either two (2 x N) or (3 x N) arrays containing pairs of half and
% difference vectors in spherical or cartesian coordinates. The function
% optionally returns the light and view directions in spherical coordinates as
% third and fourth arguments.
function [light_dirs, view_dirs, light_dirs_sph, view_dirs_sph] = rus2cart(half, diff)
    % input parsing
    if size(half, 1) == 3 && size(diff, 1) == 3
        % bring everything to spherical coordinates
        half_sph = utils.cart2sph2(half(1, :), half(2, :), half(3, :));
    elseif size(half, 1) == 2 && size(diff, 1) == 2
        half_sph = half;
        half = utils.sph2cart2(half(1, :), half(2, :));
        diff = utils.sph2cart2(diff(1, :), diff(2, :));
    elseif ~size(half, 1) == 3 && ~size(diff, 1) == 3
        error('input arrays must both be (2 x N) or (3 x N)');
    end
    n = size(half, 2);
    
    cos_theta = cos(half_sph(1, :));
    sin_theta = sin(half_sph(1, :));
    cos_phi = cos(half_sph(2, :));
    sin_phi = sin(half_sph(2, :));
    
    light_dirs = zeros(3, n);
    for jj = 1 : n
        rotZ = [cos_phi(jj), -sin_phi(jj), 0;
            sin_phi(jj), cos_phi(jj), 0;
            0, 0, 1];
        rotY = [cos_theta(jj), 0, sin_theta(jj); 
            0, 1, 0;
            -sin_theta(jj), 0, cos_theta(jj)];
        light_dirs(:, jj) = rotZ * rotY * diff(:, jj);
    end
    light_dirs = light_dirs ./ repmat(sqrt(sum(light_dirs .^ 2, 1)), 3, 1);

    % calculate wo by reflecting along the halfway vector
    view_dirs = 2. * repmat(dot(light_dirs, half), 3, 1) .* half - light_dirs;
    view_dirs = view_dirs ./ repmat(sqrt(sum(view_dirs .^ 2, 1)), 3, 1);
    
    if nargout > 2
        light_dirs_sph = utils.cart2sph2(light_dirs(1, :), light_dirs(2, :), light_dirs(3, :));
        view_dirs_sph = utils.cart2sph2(view_dirs(1, :), view_dirs(2, :), view_dirs(3, :));
    end
end
