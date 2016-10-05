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
% Given an 3 x N matrix of unit directions into the hemisphere with
% positive z-coordinate, this method adds new unit vectors that are evenly
% distributed on the unit circle in the z = 0 plane. The number of these
% added vectors depends on the distribution of the existing vectors. It is
% assumed that these vectors are aligned in planes of equal z coordinates.
% The method then tries to determine the number of vectors on the lowest
% such plane and adds the same number, or at least ten in the z = 0 plane.
function dirs = add_bottom_ring(dirs)
    assert(size(dirs, 1) == 3);
    % we're only interested in the inclination angle
    dirs_z = dirs(3, :);

    % do we actually need to add a bottom ring?
    if min(dirs_z) > 0
        % find lowest z-plane
        dirs_incl_thresh = min(dirs_z) + 1e-4 * min(dirs_z);
        dirs_lowest_plane = dirs(:, dirs_z < dirs_incl_thresh);
        num_lowest_dirs = max(10, size(dirs_lowest_plane, 2));
        
        % generate new vectors in spherical coordinates
        polar = repmat(pi / 2, 1, num_lowest_dirs);
        azimuth = linspace(0, (2 - 2 / num_lowest_dirs) * pi, num_lowest_dirs);
        dirs_bottom_ring_sph = [polar; azimuth];
        
        % convert to cartesian coordinates
        dirs = [dirs, utils.sph2cart2(dirs_bottom_ring_sph)];
    end
end
