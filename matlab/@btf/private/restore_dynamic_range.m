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
% The dynamic range of the reflectance data might have been reduced by
% storing the logarithm of the data. This facilitates matrix factorization.
function values = restore_dynamic_range(obj, values)
    if isfield(obj.meta, 'dynamic_range_reduction_method') && ...
            obj.meta.dynamic_range_reduction_method ~= 0
        
        if isfield(obj.meta, 'drr_eps')
            drr_eps = obj.meta.drr_eps;
            if isfield(obj.meta, 'drr_offset')
                drr_eps = drr_eps + obj.meta.drr_offset;
            end
        else
            drr_eps = 1e-5;
        end
        
        values = exp(values) - drr_eps;
    end
end
