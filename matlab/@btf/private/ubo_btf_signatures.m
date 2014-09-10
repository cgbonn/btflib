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
% * the terms of the GNU General Public License as published by the Free
% * Software Foundation, either version 3 of the License, or (at your
% * option) any later version.
% *
% * btflib is distributed in the hope that it will be useful, but WITHOUT
% * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% * for more details.
% *
% * You should have received a copy of the GNU General Public License along
% * with btflib.  If not, see <http://www.gnu.org/licenses/>.
% *
% *************************************************************************
%
% Generate a map that assigns a format string to each unique file signature
% string. Furthermore, for some signatures an additional flag is set that
% determines behaviour when reading the file.
function signatures_map = ubo_btf_signatures()
    % those are all currently supported file signatures for Uni Bonn BTFs
    signatures_map = containers.Map();
    signatures_map('!DFMF08FCR') =          'DFMF';
    signatures_map('!FMF12FCER') =          {'FMF', true}; % has additional field in header
    signatures_map('!FMF06FCR') =           {'FMF', false}; % those two lack
    signatures_map('!FMF06FC') =            {'FMF', false}; % the above field
    signatures_map('!PVF06FCR') =           'PVF';
    signatures_map('!PVF06FC') =            'PVF';
end