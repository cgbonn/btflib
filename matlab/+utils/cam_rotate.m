% *************************************************************************
% * Copyright 2016 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2016-01-27
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
% This function can be used to rotate an axis object's camera around the
% camera's target position. It's behavior is similar to that of Matlab's view(),
% except that it doesn't necessarily rotate around the origin. If the camera is
% panned, the center of rotation is shifted accordingly.
% Inputs are the camera's spherical coordinates. The elevation angle goes from 0
% to 180 degrees, the azimuth angle is arbitrary. Optionally, the camera's
% distance from it's point of view can be specified via the radius paramter.
% If called with only the axis object as an argument, the function only returns
% the spherical coordinates.
function [azimuth, elevation, radius] = cam_rotate(axis_handle, azimuth, elevation, radius)
    cam_target = get(axis_handle, 'CameraTarget');
    cam_pos = get(axis_handle, 'CameraPosition');
    
    % move target into origin
    cam_pos = cam_pos - cam_target;
    
    if ~exist('radius', 'var')
        radius = norm(cam_pos);
    end
    
    cam_pos_sph = utils.cart2sph2(cam_pos, true);
    
    if nargin == 1
        azimuth = rad2deg(cam_pos_sph(2));
        elevation = rad2deg(cam_pos_sph(1));
        radius = cam_pos_sph(3);
    else
        cam_pos_new = utils.sph2cart2([deg2rad(elevation); deg2rad(azimuth); radius]);
        set(axis_handle, 'CameraPosition', cam_pos_new' + cam_target);
    end
end
