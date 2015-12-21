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
% This function computes the halfway and difference vectors of the Rusinkiewicz
% parameterization. Input are either two (2 x N) or (3 x N) arrays containing
% pairs of incoming and outgoing unit length directions in spherical or
% cartesian coordinates. The function optionally returns the halfway and
% difference vectors in spherical coordinates as third and fourth arguments.
function [half, diff, half_sph, diff_sph] = cart2rus(light_dirs, view_dirs)
    % input parsing
    if size(light_dirs, 1) == 2 && size(view_dirs, 1) == 2
        % bring everything to cartesian coordinates
        light_dirs = utils.sph2cart2(light_dirs(1, :), light_dirs(2, :));
        view_dirs = utils.sph2cart2(view_dirs(1, :), view_dirs(2, :));
    elseif ~size(light_dirs, 1) == 3 && ~size(view_dirs, 1) == 3
        error('input arrays must both be (2 x N) or (3 x N)');
    end
    n = size(light_dirs, 2);
    
    % compute halfway vector & bring to spherical coordinates
    half = light_dirs + view_dirs;
    half = half ./ repmat(sqrt(sum(half .^ 2, 1)), 3, 1);
    half_sph = utils.cart2sph2(half(1, :), half(2, :), half(3, :));

    cos_theta = cos(-half_sph(2, :));
    sin_theta = sin(-half_sph(2, :));
    cos_phi = cos(-half_sph(1, :));
    sin_phi = sin(-half_sph(1, :));
    
    % compute difference vector
    diff = light_dirs;
    for jj = 1 : n,
       rot_z = [cos_theta(jj), -sin_theta(jj), 0;
           sin_theta(jj), cos_theta(jj), 0;
           0, 0, 1];
       rot_y = [cos_phi(jj), 0, sin_phi(jj);
           0, 1, 0;
           -sin_phi(jj), 0, cos_phi(jj)];
       diff(:, jj) = rot_y * rot_z * light_dirs(:, jj);
    end
    diff = diff ./ repmat(sqrt(sum(diff .^ 2, 1)), 3, 1);
    
    diff_sph = utils.cart2sph2(diff(1, :), diff(2, :), diff(3, :));
end
