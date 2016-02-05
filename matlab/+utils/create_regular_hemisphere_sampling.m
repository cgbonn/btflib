% *************************************************************************
% * Copyright 2016 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2016-01-04
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
% This function creates a regular grid in spherical coordinates on the
% hemisphere in positive z-direction. The radial component is fixed to 1.
% The first two (optional) arguments specify the number of grid points
% along the inclination and azimuth angles, and respectively default to 90
% and 360 points, i.e. appr. 1 degree grid resolution at the equator. If
% the third argument is set to true, the first and last azimuth angles will
% coincide (0 and 2 pi), which will allow to continuously interpolate
% across this boundary. The function returns the grid points in ndgrid
% format as both cartesian coordinates and optionally as spherical
% coordinates. The returned arrays respectively contain the xyz and the two
% spherical coordinates in the first dimension.
function [sampling, sampling_sph] = create_regular_hemisphere_sampling(...
        inclination_res, azimuth_res, wrap_azimuth, avoid_grazing)
    
    if ~exist('inclination_res', 'var')
        inclination_res = 90;
    end
    
    if ~exist('azimuth_res', 'var')
        azimuth_res = 360;
    end
    
    if ~exist('wrap_azimuth', 'var')
        wrap_azimuth = false;
    end
    
    if ~exist('avoid_grazing', 'var')
        avoid_grazing = false;
    end
    
    if avoid_grazing
        thetas = linspace(0, pi / 2, inclination_res + 1);
        thetas = thetas(1 : end - 1);
    else
        thetas = linspace(0, pi / 2, inclination_res);
    end
    
    if wrap_azimuth
        phis = linspace(0, 2 * pi, azimuth_res);
    else
        phis = linspace(0, 2 * pi, azimuth_res + 1);
        phis = phis(1 : end - 1);
    end
    
    [thetas, phis] = ndgrid(thetas, phis);
    
    sampling = utils.sph2cart2(thetas, phis);
    
    if nargout > 1
        sampling_sph = cat(1, reshape(thetas, [1, size(thetas)]), ...
            reshape(phis, [1, size(phis)]));
    end
end
