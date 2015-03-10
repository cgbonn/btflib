% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2015-03-10
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
% This class provides basic reading, writing and decoding functionality for
% Bidirectional Texture Functions (BTFs). The binary file formats that are
% currently supported by the class are those used in the BTFDBB at Bonn
% University (http://cg.cs.uni-bonn.de/en/projects/btfdbb/).
classdef btf < handle
    properties (GetAccess = public, SetAccess = protected)
        format_str; % string that is used to identify the BTF format
        data; % a struct that holds the components of compressed BTF files
        meta; % all other meta data stored in a struct
        
        quality; % a scalar in (0, 1] that defines the quality with which compressed btfs should be read from file
    end
    
    properties (Access = protected)
        % angular interpolation data structures
        TriL; % Delaunay triangulation of parabolic light sample coordinates
        TriV; % Delaunay triangulation of parabolic view sample coordinates
    end
    
    methods (Access = public)
        function obj = btf(file_name, quality)
            % btf constructor, which currently only supports construction by reading the data from file
            %
            % The quality flag, which is a float in (0, 1] determines how
            % much of the compressed data is read from file and thereby
            % determines the quality of reconstructed color values, e.g. in
            % the case of matrix factorizations this determines how many
            % components of a truncated SVD should be loaded.
            %
            % TODO: implement construction from data and meta data structs
            if ~exist('quality', 'var')
                obj.quality = 1;
            else
                obj.quality = quality;
            end
            
            obj = obj.read(file_name);
        end
        
        function obj = read(obj, file_name)
            % read binary BTF file
            fid = fopen(file_name, 'r');
            [obj.format_str, header_flag, signature] = identify_signature(fid);
            frewind(fid);

            switch obj.format_str
                case 'BDI'
                    [obj.data, obj.meta] = read_ubo_bdi(fid, signature, header_flag);
                case 'DFMF'
                    [obj.data, obj.meta] = read_ubo_dfmf(fid, signature, obj.quality);
                case 'FMF'
                    [obj.data, obj.meta] = read_ubo_fmf(fid, signature, header_flag, obj.quality);
                case 'PVF'
                    [obj.data, obj.meta] = read_ubo_pvf(fid, signature, obj.quality);
                otherwise
                    error(['BTF: This format is currently unsupportet, ', ...
                        'I read the signature ''%s'' and determined ', ...
                        'the format ''%s'' from it'], signature, obj.format_str);
            end

            fclose(fid);

            obj.meta.file_name = file_name;
            obj = obj.init_dirs();
        end
        
        function obj = write(obj, file_name)
            % write binary BTF file
            %
            % TODO: implement writing of PVFs
            fid = fopen(file_name, 'w');
            
            switch obj.format_str
                case 'BDI'
                    obj.write_ubo_bdi(fid);
                case 'DFMF'
                    write_ubo_dfmf(fid, obj.data, obj.meta);
                case 'FMF'
                    write_ubo_fmf(fid, obj.data, obj.meta);
                otherwise
                    error('BTF: writing is only implemented for BDI or (D)FMF format!');
            end
            
            fclose(fid);
        end
        
        function obj = write_bdi(obj, file_name)
            % write binary BDI file
            fid = fopen(file_name, 'w');
            obj.write_ubo_bdi(fid);
            fclose(fid);
        end
        
        function obj = buffer_bdi(obj, max_free_mem_percentage)
            % attempt to buffer as many ABRDF chunks from a BDI to memory;
            % optionally a fraction can be specified of the available memory
            % that is to be used
            if ~strcmp(obj.format_str, 'BDI')
                error('buffer_bdi only works with BDIs. non-BDIs are always fully loaded to memory.');
            end
            
            if ~exist('max_free_mem_percentage', 'var')
                max_free_mem_percentage = 0.9;
            end
            
            mem_available = floor(freemem() * max_free_mem_percentage);
            
            mem_required = obj.meta.chunk_size * obj.meta.num_chunks * sizeof('uint16');
            if mem_required > mem_available
                obj.data.num_chunks_in_buffer = floor(obj.meta.num_chunks * mem_available / mem_required);
            else
                obj.data.num_chunks_in_buffer = obj.meta.num_chunks;
            end
            
            if mem_available < obj.meta.chunk_size * sizeof('uint16')
                obj.data.num_chunks_in_buffer = 0;
                error('too little memory available for buffering, one chunk requires %3.2f MB', obj.meta.chunk_size * sizeof('uint16') / 1024 / 1024);
            end
            
            % allocate memory
            obj.data.chunks = zeros(obj.meta.chunk_size, ...
                obj.data.num_chunks_in_buffer, 'uint16');
            
            % read chunks from file
            obj.data.fid = fopen(obj.meta.file_name, 'r');
            for c = 1 : obj.data.num_chunks_in_buffer - 1
                fprintf('reading chunk %d / %d (total: %d)\n', ...
                    c, obj.data.num_chunks_in_buffer, obj.meta.num_chunks);
                obj.get_bdi_chunk(c);
            end
            c = obj.data.num_chunks_in_buffer;
            % handle last chunk separately because it is usually smaller, but
            % indexing the memory chunks in each of the above iterations would
            % strongly reduce performance
            obj.get_bdi_chunk(c);
            fclose(obj.data.fid);
            obj.data.chunks_buffered(1 : obj.data.num_chunks_in_buffer) = true;
        end
        
        function obj = only_use_buffered(obj, value)
            % if a bdi is partially buffered, setting this to true will tell the
            % decoder to ignore those pixels that aren't buffered, otherwise, it
            % will (slowly) read those missing pixel from file
            if ~strcmp(obj.format_str, 'BDI')
                error('this only works with BDIs. non-BDIs are always fully loaded to memory.');
            end
            obj.data.only_use_buffered = logical(value(1));
        end
        
        function obj = set_data(obj, data)
            % very simple setter method for changing the compressed BTF data
            % TODO: perform sanity checks depending on the format!
            obj.data = data;
        end
        
        function obj = set_meta(obj, varargin)
            % setter method for changing meta data
            %
            % arguments can either be a single struct containing all
            % necessary fields, or an arbitrary number of key-value pairs
            % that will be used to update or add fields in the existing
            % meta data struct
            if numel(varargin) == 1
                % replace the full meta data struct
                obj.meta = varargin{1};
            else
                % we expect key-value pairs in varargin
                assert(mod(numel(varargin), 2) == 0);
                for a = 1 : 2 : numel(varargin)
                    key = varargin{a};
                    value = varargin{a + 1};
                    obj.meta.(key) = value;
                end
            end
        end
        
        function [L, V] = get_LV_full(obj)
            % get cartesian product of the sampled light and view directions
            L = obj.get_L_full();
            V = obj.get_V_full();
        end
        
        function L = get_L_full(obj)
            % get light direction for each sample in the data
            L = repmat(obj.meta.L, obj.meta.nV, 1);
        end
        
        function V = get_V_full(obj)
            % get view direction for each sample in the data
            V = obj.meta.V(:) * ones(1, obj.meta.nL);
            V = reshape(V', [], 3);
        end
        
        function [lin, laz, vin, vaz] = inds_to_angles(obj, l, v)
            % given light and view indices, return the light inclination and
            % azimuth, as well as the view inclination and azimuth angles
            l_sph = cart2sph2(obj.meta.L(l, :));
            v_sph = cart2sph2(obj.meta.V(v, :));
            lin = l_sph(1);
            laz = l_sph(2);
            vin = v_sph(1);
            vaz = v_sph(2);
        end
        
        function abrdf = decode_abrdf(obj, x, y)
            % decode a full ABRDF for a given pair of x- and y-coordinates
            %
            % no interpolation yet, so x and y should be integers
            x = uint32(x);
            y = uint32(y);
            switch obj.format_str
                case 'BDI'
                    abrdf = obj.decode_bdi_abrdf(x, y);
                case 'DFMF'
                    abrdf = obj.decode_dfmf_abrdf(x, y);
                case 'FMF'
                    abrdf = obj.decode_fmf_abrdf(x, y);
                case 'PVF'
                    abrdf = obj.decode_pvf_abrdf(x, y);
                otherwise
                    error('decoder for BTF in %s format not implemented!', obj.format_str);
            end
            
            % restore RGB values from potentially decorrelated data
            abrdf = obj.undecorrelate(abrdf);
            
            % invert dynamic range reduction if necessary
            abrdf = obj.restore_dynamic_range(abrdf);
            
            % remove cosine term, if it is in the data
            if obj.meta.cosine_flag
                abrdf = abrdf ./ repmat(obj.meta.L(:, 3), 1, obj.meta.nV, obj.meta.num_channels);
            end
            
            % rearrange as image
            abrdf = reshape(abrdf, [obj.meta.nL, obj.meta.nV, obj.meta.num_channels]);
        end
        
        function texel = decode_texel(obj, x, y, L, V)
            % decode a single texel at given pixel coordinates for given light & view directions
            % 
            % x, y need to be integer indices into the texture space, L and
            % V can either be integers indices in the angular domain or two
            % vectors (either 3D cartesian or 2D inclination-azimuth pairs).
            % In the latter case, angular interpolation is performed
            % between the color values corresponding to the closest sampled
            % direction pairs.
            if isscalar(L) && isscalar(V)
                L = uint32(L);
                V = uint32(V);
                % direct access via indices
                switch obj.format_str
                    case 'BDI'
                        texel = obj.decode_bdi_texel(x, y, L, V);
                    case 'DFMF'
                        texel = obj.decode_dfmf_texel(x, y, L, V);
                    case 'FMF'
                        texel = obj.decode_fmf_texel(x, y, L, V);
                    case 'PVF'
                        texel = obj.decode_pvf_texel(x, y, L, V);
                    otherwise
                        error('decoder for BTF in %s format not implemented!', obj.format_str);
                end
                
                % restore RGB values from potentially decorrelated data
                texel = obj.undecorrelate(texel);

                % invert dynamic range reduction if necessary
                texel = obj.restore_dynamic_range(texel);

                % remove cosine term, if it is in the data
                if obj.meta.cosine_flag
                    texel = texel ./ obj.meta.L(L, 3);
                end
            elseif numel(L) == 2 && numel(V) == 2 || numel(L) == 3 && numel(V) == 3
                % access via light and view directions -> interpolate samples
                texel = zeros(obj.nC, 1, obj.data.class);
                
                % get interpolation indices and weights
                [ls, vs, weights_l, weights_v] = obj.lookup_dirs(L, V);
                for li = 1 : numel(ls)
                    l = ls(li);
                    for vi = 1 : numel(vs)
                        v = vs(vi);
                        
                        texel = texel + weights_l(li) * weights_v(vi) * ...
                            obj.decode_texel(x, y, l, v);
                    end
                end
            else
                error('error while decoding: either provide readily computed indices, 3D cartesian directions, or two polar angles respectively for light and view.');
            end
        end
        
        function img = decode_texture(obj, L, V)
            % decode full texture for given light & view directions
            %
            % L and V can either be integers indices in the angular domain
            % or two vectors (either 3D cartesian or 2D inclination-azimuth
            % pairs). In the latter case, angular interpolation is
            % performed between the color values corresponding to the
            % closest sampled direction pairs.
            if isscalar(L) && isscalar(V)
                L = uint32(L);
                V = uint32(V);
                % direct access via indices
                switch obj.format_str
                    case 'BDI'
                        img = obj.decode_bdi_texture(L, V);
                    case 'DFMF'
                        img = obj.decode_dfmf_texture(L, V);
                    case 'FMF'
                        img = obj.decode_fmf_texture(L, V);
                    case 'PVF'
                        img = obj.decode_pvf_texture(L, V);
                    otherwise
                        error('decoder for BTF in %s format not implemented!', obj.format_str);
                end
                
                % restore RGB values from potentially decorrelated data
                img = obj.undecorrelate(img);

                % invert dynamic range reduction if necessary
                img = obj.restore_dynamic_range(img);

                % remove cosine term, if it is in the data
                if obj.meta.cosine_flag
                    img = img ./ obj.meta.L(L, 3);
                end
                
                % rearrange as image
                img = reshape(img, [obj.meta.width, obj.meta.height, obj.meta.num_channels]);
            elseif numel(L) == 2 && numel(V) == 2 || numel(L) == 3 && numel(V) == 3
                % access via light and view directions -> interpolate samples
                img = zeros(obj.meta.width, obj.meta.height, obj.meta.num_channels, obj.data.class);
                [ls, vs, weights_l, weights_v] = obj.lookup_dirs(L, V);
                for li = 1 : numel(ls)
                    l = ls(li);
                    for vi = 1 : numel(vs)
                        v = vs(vi);
                        
                        img = img + weights_l(li) * weights_v(vi) * ...
                            obj.decode_texture(l, v);
                    end
                end
            else
                error('error while decoding: either provide readily computed indices, 3D cartesian directions, or two polar angles respectively for light and view.');
            end
        end
        
        function tensor = get_6d_tensor(obj)
            % TODO: implement extraction of the full data tensor arranged as a
            % 6-dimensional array (or 7D if the color channels aren't unrolled)
            error('not implemented yet!');
        end
    end
    
    methods (Access = protected)
        function obj = init_dirs(obj)
            % set up data structures for interpolation in angular domain
            % (Delaunay triangulations in parabolic coordinates) for
            % performing angular lookup for arbitrary light / view
            % directions
            obj.meta.L = double(obj.meta.L);
            obj.meta.V = double(obj.meta.V);
            L = obj.meta.L;
            V = obj.meta.V;
            
            % we potentially have to add a bottom ring in case our angular
            % sampling doesn't reach low enough in the hemisphere (we'd
            % then have a margin around our triangulation in which we
            % couldn't determine any enclosing triangles for queried
            % directions)
            L = add_bottom_ring(L);
            V = add_bottom_ring(V);
            
            % set up Delaunay trianglulation of parabolic map of both light
            % and view directions
            L_par = cart2par(L);
            V_par = cart2par(V);
            obj.TriL = delaunayTriangulation(L_par);
            obj.TriV = delaunayTriangulation(V_par);
        end
        
        function [L_idxs, V_idxs, L_weights, V_weights] = lookup_dirs(obj, L, V)
            % determine indices and weights for angular interpolation
            
            % transform input directions to parabolic coordinates
            if numel(L) == 3 && numel(V) == 3
                L_par = cart2par(L);
                V_par = cart2par(V);
            elseif numel(L) == 2 && numel(V) == 2
                % we assume L and V to be in parabolic coordinates already!
                % unfortunately, there is no sane way to distinguish those
                % from spherical coordinates
                L_par = L;
                V_par = V;
            end
            
            % lookup in Delaunay triangulation for both directions
            [L_ti, L_bc] = obj.TriL.pointLocation(L_par);
            [V_ti, V_bc] = obj.TriV.pointLocation(V_par);
            if any(isnan(L_ti)) || any(isnan(V_ti))
                % one of the queried directions is outside the
                % triangulation, i.e. below the surface -> we don't need to
                % interpolate at all
                L_idxs = [];
                V_idxs = [];
                L_weights = [];
                V_weights = [];
                return;
            end
            
            % we cannot return indices into the bottom ring, as we don't
            % have any data available there
            L_idxs_with_bottom_ring = obj.TriL.ConnectivityList(L_ti, :);
            V_idxs_with_bottom_ring = obj.TriV.ConnectivityList(V_ti, :);
            
            % hence, we're only interested in those indices and weights
            % that are actually containted in the sampling
            L_idxs = L_idxs_with_bottom_ring(L_idxs_with_bottom_ring < obj.meta.nL);
            V_idxs = V_idxs_with_bottom_ring(V_idxs_with_bottom_ring < obj.meta.nV);
            L_weights = L_bc(L_idxs_with_bottom_ring < obj.meta.nL);
            V_weights = V_bc(V_idxs_with_bottom_ring < obj.meta.nV);
        end
    end
end
