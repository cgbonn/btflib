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
% hemisphere in positive z-direction. The radial component is fixed to 1. The
% two (optional) arguments specify the number of grid points along the
% inclination and azimuth angles, and respectively default to 90 and 180 points,
% i.e. appr. 1 degree grid resolution at the equator. The function returns the
% grid points in ndgrid format as both cartesian coordinates and optionally as
% spherical coordinates. The returned arrays respectively contain the xyz and
% the two spherical coordinates in the first dimension.
%
function [sampling, sampling_sph] = create_regular_hemisphere_sampling(...
        inclination_res, azimuth_res)
    
    if ~exist('inclination_res', 'var')
        inclination_res = 90;
    end
    
    if ~exist('azimuth_res', 'var')
        azimuth_res = 180;
    end
    
    thetas = linspace(0, pi / 2, inclination_res);
    phis = linspace(0, (2 - 2 / azimuth_res) * pi, azimuth_res);
    
    [thetas, phis] = ndgrid(thetas, phis);
    
    sampling = utils.sph2cart2(thetas, phis);
    
    if nargout > 1
        sampling_sph = cat(1, reshape(thetas, [1, size(thetas)]), ...
            reshape(phis, [1, size(phis)]));
    end
end
