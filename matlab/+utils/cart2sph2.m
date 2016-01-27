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
% Replacement for Matlab's cart2sph with the following differences:
%   - instead of an elevation, this function computes a polar angle, i.e.
%     the angular distance in radians from the zenith (x = 0, y = 0, z = 1)
%   - per default, the radial coordinate is not computed as this function
%     is intended for normalized vectors; if it is required, an additional
%     boolean input parameter needs to be set to true
%   - input can either be a single N x 3 or 3 x N matrix which stores the
%     cartesian x-, y- and z-coordinates in its rows or columns, or separate N x
%     1 or 1 x N arrays of these three quantities
%   - by default it either returns an N x 2 or 2 x N array with the polar angle
%     in its first and the azimuth angle in its second column / row; or if
%     requested, it returns an N x 3 or 3 x N array which furthermore stores the
%     radial component in its third column
%
% Usage:
%   sph = cart2sph2(xyz)
%   sph = sph2cart2(xyz, with_radius)
%   sph = sph2cart2(x, y, z)
%   sph = sph2cart2(x, y, z, with_radius)
% 
% where all input arrays are N x 3, 3 x N, N x 1 or 1 x N respectively, the
% optional last argument is a boolean and the output array is N x 2 or 2 x N.
function sph = cart2sph2(varargin)
    % input parsing
    if numel(varargin) < 3
        % one input matrix of dimensions n x 3
        if size(varargin{1}, 1) == 3
            x = utils.sel_dim(varargin{1}, 1, 1);
            y = utils.sel_dim(varargin{1}, 1, 2);
            z = utils.sel_dim(varargin{1}, 1, 3);
            n = size(x);
            n(1) = [];
        else
            x = utils.sel_dim(varargin{1}, 2, 1);
            y = utils.sel_dim(varargin{1}, 2, 2);
            z = utils.sel_dim(varargin{1}, 2, 3);
            n = size(x);
            n(2) = [];
        end
    elseif numel(varargin) >= 3
        % three separate matrices
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        n = size(x);
        if numel(n) > 1 && n(1) == 1
            n = n(2 : end);
        end
    end
    
    if numel(varargin) > 3
        with_radius = varargin{4};
    else
        if numel(varargin) == 2
            with_radius = varargin{2};
        else
            with_radius = false;
        end
    end
    
    len_xy = hypot(x, y);
    
    % radius?
    if with_radius
        sph = zeros([3, prod(n)], class(x)); %#ok<ZEROLIKE>
        
        sph(3, :) = reshape(hypot(len_xy, z), 1, []);
    else
        sph = zeros([2, n], class(x)); %#ok<ZEROLIKE>
    end
    
    % inclination
    sph(1, :) = reshape(atan2(len_xy, z), 1, []);
    
    % azimuth
    sph(2, :) = reshape(atan2(y, x), 1, []);
    
    % let's try to maintain the array orientation
    if numel(n) == 2 && any(n == 1) || numel(n) == 1
        % special case: 2D array -> keep orientation
        if numel(varargin) < 3 && (size(varargin{1}, 2) == 3 && size(varargin{1}, 1) ~= 3) || ...
                numel(varargin) >= 3 && size(varargin{1}, 2) == 1 && ~all(n == 1)
            sph = sph';
        else
            if numel(n) ~= 1
                sph = reshape(sph, [size(sph, 1), n(2)]);
            end
        end
    else
        % otherwise we default to putting the xyz-coordinates into the first
        % dimension and keeping the the shape of the input
        sph = reshape(sph, [size(sph, 1), n]);
    end
end
