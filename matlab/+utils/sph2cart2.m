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
%   - input can either be a single N x 3 or 3 x N or N x 2 or 2 x N matrix which
%     stores the polar angle in its first, the azimuth angle in its second and
%     optionally the radius in its third column or row, or separate N x 1 or 1 x
%     N arrays of these three quantities
%   - it returns an N x 3 or 3 x N array of cartesian coordinates
%
% Usage:
%   xyz = sph2cart2(sph)
%   xyz = sph2cart2(polar, azimuth)
%   xyz = sph2cart2(polar, azimuth, radius)
% 
% where all input and out arrays are N x 3, 3 x N, N x 1 or 1 x N respectively.
function xyz = sph2cart2(varargin)
    % input parsing
    if numel(varargin) == 1
        if size(varargin{1}, 2) == 2 || size(varargin{1}, 2) == 3
            % one input matrix of dimensions N x 2 or N x 3
            inclination = varargin{1}(:, 1);
            azimuth = varargin{1}(:, 2);
            n = size(inclination);

            if size(varargin{1}, 2) == 3
                radius = varargin{1}(:, 3);
            else
                radius = ones(n, class(inclination)); %#ok<ZEROLIKE>
            end
        else
            % one input matrix of dimensions 2 x N or 3 x N
            inclination = varargin{1}(1, :);
            azimuth = varargin{1}(2, :);
            n = size(inclination);
            
            if size(varargin{1}, 1) == 3
                radius = varargin{1}(3, :);
            else
                radius = ones(n, class(inclination)); %#ok<ZEROLIKE>
            end
        end
    elseif numel(varargin) >= 2
        % separate input vectors
        inclination = varargin{1};
        azimuth = varargin{2};
        n = size(inclination);
        
        if numel(varargin) >= 3
            radius = varargin{3};
        else
            radius = ones(size(inclination), class(inclination)); %#ok<ZEROLIKE>
        end
    end
    
    xyz = zeros([3, prod(n)], class(azimuth)); %#ok<ZEROLIKE>
    
    rsinelev = radius .* sin(inclination);
    xyz(1, :) = reshape(rsinelev .* cos(azimuth), 1, []);
    xyz(2, :) = reshape(rsinelev .* sin(azimuth), 1, []);
    xyz(3, :) = reshape(radius .* cos(inclination), 1, []);
    
%     % let's try to maintain the array orientation
%     di = find(n == 1, 1); % check if inputs were vectors or one 2D array
%     if numel(n) == 2 && ~isempty(di)
%         % special case: 2D array -> keep orientation
%         if di == 2
%             xyz = xyz';
%         end
%     else
%         % otherwise we default to putting the xyz-coordinates into the first
%         % dimension and keeping the the shape of the input
%         xyz = reshape(xyz, [3, n]);
%     end
    
    % let's try to maintain the array orientation
    if numel(n) == 2 && any(n == 1)
        % special case: 2D array -> keep orientation
        if (numel(varargin) == 1 && size(varargin{1}, 2) == 2 || ...
                size(varargin{1}, 2) == 3 || ...
                numel(varargin) >= 2 && size(varargin{1}, 2) == 1)
            xyz = xyz';
        else
            xyz = reshape(xyz, [3, n(2)]);
        end
    else
        % otherwise we default to putting the xyz-coordinates into the first
        % dimension and keeping the the shape of the input
        xyz = reshape(xyz, [3, n]);
    end
end
