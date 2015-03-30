% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-03-30
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
% 
function [L, V] = create_ubo_dome1_sampling()
    inclination_angles_degree = [0, 11, 23.5, 30, 37.5, 45, 52.5, 60, 67.5, 75];
    azimuth_ring_count = [1, 6, 12, 12, 12, 18, 18, 24, 24, 24];
    azimuth_ring_offsets_degree = [0, 15, 0, 15, 0, 0, 7.5, 0, 7.5, 0];
    
    n_dirs = sum(azimuth_ring_count);
    inclinations = zeros(n_dirs, 1);
    azimuths = zeros(n_dirs, 1);
    % top dir
    inclinations(1) = inclination_angles_degree(1);
    azimuths(1) = azimuth_ring_offsets_degree(1);
    
    % iterate over rings
    ris = cumsum(azimuth_ring_count);
    for r = 2 : numel(azimuth_ring_count)
        ri = ris(r - 1) + 1 : ris(r);
        inclinations(ri) = inclination_angles_degree(r);
        
        azims = linspace(0, 360 * (1 - 1 / azimuth_ring_count(r)), azimuth_ring_count(r));
        azims = azims + azimuth_ring_offsets_degree(r);
        azims(azims >= 180) = mod(azims(azims >= 180), 180) - 180;
        azimuths(ri) = azims;
    end
    
    L = utils.sph2cart2([deg2rad(inclinations), deg2rad(azimuths)]);
    V = L;
end
