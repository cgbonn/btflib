% *************************************************************************
% * Copyright 2016 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2016-01-09
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
% Normalize a set of vectors according to their Euclidean norm. The user
% can either explicitly specify a dimension along which the norms should
% be computed, or it can be guessed automatically (in this case, the
% function works along the smaller dimension, except for the case of a
% single vector).
function vectors = normalize(vectors, dim)
    tmp = size(vectors);
    
    % try to automatically determine along which dimension we want to normalize
    if ~exist('dim', 'var')
        if min(tmp) == 1
            % we have a single vector
            [~, dim] = max(tmp);
        else
            % special case: we have a set of 2D or 3D vectors
            [m, mi] = min(tmp);
            if m == 2 || m == 3
                [~, dim] = setdiff(tmp, mi);
                dim = dim(1);
            else
                % by default, work along the smaller dimension
                dim = mi;
            end
        end
    end

    if dim < 1 || dim > 2 || numel(tmp) > 2
        error( 'only vectors or matrices are supported' );
    end
    
    tmp(setdiff(1 : 2, dim)) = 1;
    
    vectors = vectors ./ repmat(sqrt(sum(vectors .^ 2, dim)), tmp);
end
