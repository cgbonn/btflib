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
% This function converts parabolic coordinates specified as an N x 2 or 2 x N
% array into 3D cartesian coordinates according to the relations l_{par} =
% x_{par} ^ 2 + y_{par}; xy_{car} = 2 * xy_{par} / (l + 1); z_{car} = (1 - l) /
% (l + 1) and returns the cartesian x-, y- and z-coordinates in an N x 3 or 3 x
% N array.
%
% Usage:
%   xyz_car = par2cart(xy_par)
%   xyz_car = par2cart(x_par, y_par, z_par)
% 
% where all input arrays are N x 2, 2 x N, N x 1 or 1 x N respectively and the
% output array is N x 3 or 3 x N.
function car = par2cart(varargin)
    % input parsing
    if numel(varargin) < 2
        % one input matrix of dimensions n x 3
        if size(varargin{1}, 1) == 2
            x = varargin{1}(1, :);
            y = varargin{1}(2, :);
        else
            x = varargin{1}(:, 1);
            y = varargin{1}(:, 2);
        end
    elseif numel(varargin) >= 2
        % three separate matrices
        x = varargin{1};
        y = varargin{2};
    end
    n = size(x);
    
    car = zeros([3, n], class(x)); %#ok<ZEROLIKE>
    
    length_sqr = x .^ 2 + y .^ 2;
    car(1, :) = reshape(2 * x ./ (length_sqr + 1), 1, []);
    car(2, :) = reshape(2 * y ./ (length_sqr + 1), 1, []);
    car(3, :) = reshape((1 - length_sqr) ./ (length_sqr + 1), 1, []);
    
    % let's try to maintain the array orientation
    if numel(n) == 2 && any(n == 1)
        % special case: 2D array -> keep orientation
        if (numel(varargin) == 1 && size(varargin{1}, 2) == 2 || ...
                numel(varargin) >= 2 && size(varargin{1}, 2) == 1)
            car = car';
        else
            car = reshape(car, [3, n(2)]);
        end
    else
        % otherwise we default to putting the xyz-coordinates into the first
        % dimension and keeping the the shape of the input
        car = reshape(car, [3, n]);
    end
end
