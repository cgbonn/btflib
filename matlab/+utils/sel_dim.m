% *************************************************************************
% * Copyright 2016 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2016-01-06
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
% This function can be used to index any dimension of an arbitrarily shaped
% array, without changing the size of the returned array.
%
% Usage: values = sel_dim(values, dim, index), where
%
% - values is an arbitrarily shaped array
% - dim is an integer specifying the dimension which should be indexed
% - index is an array of indices that are used to index dimension dim of the
%   input array
%
% returns
% - an array with the same shape as the input array except for the dimension
%   that was indexed
function values = sel_dim(values, dim, index)
    nd = ndims(values);
    
    index_set = repmat({':'}, 1, nd);
    index_set{dim} = index;
    
    values = values(index_set{:});
end
