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
% Read signature string from file and use it to determine the file format.
% For this, a map of supported signature strings is generated that assigns
% a format string to each signature.
function [format_str, header_flag, signature] = identify_signature(fid)
    signatures_map = ubo_btf_signatures();
    
    % read signature
    signatures = signatures_map.keys;
    max_sig_len = max(cellfun(@length, signatures));
    sig_str = fread(fid, max_sig_len, '*char');
    
    % iterate through all signatures and compare with the one in the file
    format_str = '';
    max_len_matched = 0;
    for s=1:numel(signatures)
        sig = signatures{s};
        if strcmp(sig, sig_str(1:numel(sig))')
            % we want the _longest_ matching string!
            if numel(sig) > max_len_matched
                signature = sig;
                format_str = signatures_map(sig);
                header_flag = false;
                
                % some signature strings have additionally assigned a flag
                % that further determines behaviour when reading the file
                if iscell(format_str)
                    header_flag = format_str{2};
                    format_str = format_str{1};
                end
            end
        end
    end
    
    if isempty(format_str)
        error('error determining the BTF format, I read the following signature: %s', sig_str);
    end
end