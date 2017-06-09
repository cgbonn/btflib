% *************************************************************************
% * Copyright 2015 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2016-12-30
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
% This class also allows to create BTF objects from data already stored in
% memory.
classdef btf < handle
    properties (GetAccess = public, SetAccess = protected)
        format_str; % string that is used to identify the BTF format
        data; % a struct that holds the components of compressed BTF files
        meta; % all other meta data stored in a struct
        verbose = true; % enables some status output, when lengthy operations are performed
        
        quality; % a scalar in (0, 1] that defines the quality with which compressed btfs should be read from file
    end
    
    properties (Access = protected)
        TriL; % Delaunay triangulation of parabolic light sample coordinates used for angular interpolation
        TriV; % Delaunay triangulation of parabolic view sample coordinates used for angular interpolation
        
        progress_fcn; % function handle to display updates for lengthy processes
    end
    
    methods (Access = public)
        function obj = btf(varargin)
            % btf constructor, which either reads data from file or directly from memory
            %
            % obj = btf(file_name, ...)
            % obj = btf(format_str, meta_struct, data_struct, ...)
            %
            % When reading from file, the first argument is a file name or a
            % full path to a file.
            %
            % When the object is constructed from data in memory, the first
            % argument is a string specifying the BTF format. Currently
            % supported are BDI, DFMF, FMF, PVF. The following two arguments are
            % structs containing first the meta data and second the actual BTF
            % data, which is either provided as a full BDI tensor (or
            % re-arranged matrix), or as the components of a compressed BTF
            % matrix (FMF, DFMF, PVF), e.g. the U and SxV matrices of a FMF. The
            % exact required format of  the structs is specified below.
            %
            % There are some general parameters that can be specified as
            % name-value pairs:
            %
            % The 'quality' parameter, which is a float in (0, 1] determines how
            % much of the compressed data is read from file and thereby
            % determines the quality of reconstructed color values, e.g. in the
            % case of matrix factorizations this determines how many components
            % of a truncated SVD should be loaded.
            % It is possible to provide a function handle that is called to
            % provide progress information to the user for lengthy processes.
            % This can be achieved by providing a function handle for the
            % 'progress_callback' parameter.
            %
            % Format of the meta data struct, listed are the field names, their
            % default values in braces, the requirements for their values and
            % in square brackets the BTF formats the fields are applicable to
            % (all, if nothing is specified):
            %
            % width (-1): isscalar(x) && isnumeric(x) && x > 0
            % height (-1): isscalar(x) && isnumeric(x) && x > 0
            % num_lights (151): isscalar(x) && isnumeric(x) && x > 0
            % num_views (151): isscalar(x) && isnumeric(x) && x > 0
            % num_channels (3), isscalar(x) && isnumeric(x) && x > 0
            % light_dirs ([]): isnumeric(x) && size(x, 2) == 3
            % view_dirs ([]): isnumeric(x) && size(x, 2) == 3
            % cosine_flag (false): isscalar(x) && islogical(x)
            % dynamic_range_reduction_method (0): isscalar(x) && isnumeric(x)
            % color_model (0): isnumeric(x) && isscalar(x); [DFMF]; see btf.undecorrelate
            % color_mean ([0, 0, 0]): isnumeric(x); [DFMF]
            % color_transformation_matrix (eye(3)): isnumeric(x); [DFMF]
            %
            % If num_lights == 151 and num_views == 151, light_dirs and
            % view_dirs become optional, as the default UBO Dome 1 sampling is
            % applied if no light or view directions are provided.
            %
            % The format of the data struct depends on the BTF format:
            % BDI:
            % chunks (required): this array either needs to be a tensor with
            %	dimensions nC x nL x nV x w x h or a matrix with dimensions (nC
            %	* nL * nV) x (w * h) or a matrix with dimensions (nC * nL * nV *
            %	w * S) x P, where S is the number of scan lines per data chunk
            %	and P equals ceil(h / S)
            %
            % DFMF:
            % U (required): cell array with number of elements equal to the
            %	channels in the used color model; each element of U must be a (h
            %	* w) x P or a (nL * nV) x P matrix, where P is the number of PCA
            %	components in each channel
            % SxV (required): cell array as above, only that the matrix
            %	dimensions for each element are the respectively other than
            %	those in U (if U is (h * w) x P, than SxV must be (nL * nV) x P)
            % S (optional): P x 1 array containing the singular values (unused
            %	at the moment and can thus be savely omitted)
            %
            % FMF:
            % U (required): matrix with dimensions (nC * w * h) x P or (nC * nL
            %	* nV) x P, where P is the number of PCA components
            % SxV (required): matrix with dimensions: if U is (nC * w * h) x P,
            %	then SxV must be (w * h) x P, if U is (nC * nL * nV) x P, then
            %	SxV must be (w * h) x P
            % S (optional): P x 1 array of singular values(unused at the moment
            %	and can thus be savely omitted)
            %
            % PVF:
            % Cs (required): cell array with nV elements, where each element is
            %	a (nC * nL) x P or a (nC * h * w) x P matrix, where P is the
            %	number of PCA components per view slot; those are the PCA
            %	coefficients
            % Ws (required): cell array with nV elements, where each element is
            %	a matrix; if Cs{v} is a (nC * h * w) x P matrix, then Ws{v} must
            %	be a nL x P matrix; if Cs{v} is a (nC * nL) x P matrix, then
            %	Ws{v} must be a (w * h) x P matrix; those are the PCA weights
            % Ms (required): cell array with nV elements, where each element is
            %	a vector whose length corresponds to the number of rows of the
            %	elements of Cs (if e.g. Cs{v} is a (nC * nL) x P matrix, then
            %	Ms{v} must be a vector with nC * nL elements); those are the
            %	column means of each view slot matrix
            % EVs (optional): cell array with vectors containing the eigenvalues
            %	for each view slot
            p = inputParser;
            p.CaseSensitive = false;
            p.KeepUnmatched = true;
            p.StructExpand = false;
            p.addParameter('quality', 1, @(x) isscalar(x) && isnumeric(x));
            p.addParameter('progress_callback', @obj.default_progress_fcn, ...
                @(x) isa(x, 'function_handle'));
            
            [~, supported_formats] = ubo_btf_signatures();
            
            file_name = varargin{1};
            if ~ischar(file_name)
                error('The first input argument must be a string (file name or format string).');
            end
            if ~ismember(lower(file_name), lower(supported_formats))
                % supposedly the path to a BTF file
                if exist(file_name, 'file') == 7
                    % input is a folder
                    obj.format_str = 'images';
                elseif exist(file_name, 'file') ~= 2
                    % input is neither file, nor folder
                    error('BTF file not found.');
                end
                varargin = varargin(2 : end);
                % convert raw image format to BDI?
                p.addParameter('to_bdi', true, @isscalar);
                p.parse(varargin{:});
                obj.quality = p.Results.quality;
                obj.progress_fcn = p.Results.progress_callback;
                
                % read from file
                obj = obj.read(file_name, p.Unmatched);
                
                if strcmpi(obj.format_str, 'images') && p.Results.to_bdi
                    % currently the only option to decode from raw images
                    obj.format_str = 'bdi';
                    obj.create_ubo_bdi(obj.meta, obj.data);
                end
            elseif ismember(lower(file_name), lower(supported_formats))
                % create btf object directly from data
                
                p.addRequired('format_str', @(x) any(validatestring(x, supported_formats)));
                p.addRequired('meta', @isstruct);
                p.addRequired('data', @(x) isstruct(x) || isnumeric(x));
                p.parse(varargin{:});
                obj.format_str = p.Results.format_str;
                obj.quality = p.Results.quality;
                obj.progress_fcn = p.Results.progress_callback;
                
                obj.meta.file_name = '';
                
                switch lower(obj.format_str)
                    case 'bdi'
                        obj.create_ubo_bdi(p.Results.meta, p.Results.data);
                    case 'dfmf'
                        obj.create_ubo_dfmf(p.Results.meta, p.Results.data);
                    case 'fmf'
                        obj.create_ubo_fmf(p.Results.meta, p.Results.data);
                    case 'pvf'
                        obj.create_ubo_pvf(p.Results.meta, p.Results.data);
                    otherwise
                        error('format_str %s not supported', obj.format_str);
                end
                
                obj = obj.init_dirs();
            end
        end
        
        function obj = read(obj, file_name, varargin)
            % load binary BTF file
            if ~strcmpi(obj.format_str, 'images')
                fid = fopen(file_name, 'r');
                [obj.format_str, header_flag, signature] = identify_signature(fid);
                frewind(fid);
            end

            switch lower(obj.format_str)
                case 'images'
                    [obj.data, obj.meta] = read_ubo_images(file_name, varargin{:});
                case 'bdi'
                    [obj.data, obj.meta] = read_ubo_bdi(fid, signature, header_flag);
                case 'dfmf'
                    [obj.data, obj.meta] = read_ubo_dfmf(fid, signature, obj.quality);
                case 'fmf'
                    [obj.data, obj.meta] = read_ubo_fmf(fid, signature, header_flag, obj.quality);
                case 'pvf'
                    [obj.data, obj.meta] = read_ubo_pvf(fid, signature, obj.quality);
                otherwise
                    error(['BTF: This format is currently unsupportet, ', ...
                        'I read the signature ''%s'' and determined ', ...
                        'the format ''%s'' from it'], signature, obj.format_str);
            end

            if ~strcmpi(obj.format_str, 'images')
                fclose(fid);
            end

            % store full path
            fn = which(file_name);
            if isempty(fn)
                % path already is absolute and not on the search path
                obj.meta.file_name = file_name;
            else
                obj.meta.file_name = which(file_name);
            end
            obj = obj.init_dirs();
        end
        
        function obj = write(obj, file_name)
            % save binary BTF file
            %
            % TODO: implement writing of PVFs
            fid = fopen(file_name, 'w');
            
            switch lower(obj.format_str)
                case 'bdi'
                    obj.write_ubo_bdi(fid);
                case 'dfmf'
                    write_ubo_dfmf(fid, obj.data, obj.meta);
                case 'fmf'
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
        
        function obj = crop(obj, roi, strides, output_file_name)
            % apply a region of interest to a BTF or BDI
            if ~exist('strides', 'var')
                strides = [1; 1];
            end
            
            roi = round(roi);
            strides = round(strides);
            
            % assert roi is compatible with BTF dimensions
            roi(1, 1) = max(1, min(obj.meta.width, roi(1, 1)));
            roi(1, 2) = max(roi(1, 1), min(obj.meta.width, roi(1, 2)));
            roi(2, 1) = max(1, min(obj.meta.height, roi(2, 1)));
            roi(2, 2) = max(roi(2, 1), min(obj.meta.height, roi(2, 2)));
            
            % and also ensure that the strides are not too large
            roi_dims = roi(:, 2) - roi(:, 1) + 1;
            strides(1) = max(1, min(roi_dims(1), strides(1)));
            strides(2) = max(1, min(roi_dims(2), strides(2)));
            inds_x = roi(1, 1) : strides(1) : roi(1, 2);
            inds_y = roi(2, 1) : strides(2) : roi(2, 2);
            
            % apply roi to data matrices
            switch lower(obj.format_str)
                case 'bdi'
                    obj = obj.crop_bdi(roi, strides, output_file_name);
                case 'dfmf'
                    obj = obj.crop_dfmf(roi, strides);
                case 'fmf'
                    obj = obj.crop_fmf(roi, strides);
                case 'pvf'
                    obj = obj.crop_pvf(roi, strides);
            end
            
            obj.meta.width = numel(inds_x);
            obj.meta.height = numel(inds_y);
        end
        
        function bdi = is_bdi(obj)
            % check if the object is a uncompressed BDI
            bdi = strcmpi(obj.format_str, 'bdi');
        end
        
        function tf = is_spectral(obj)
            % returns true if the BTF object stores a wavelength sampling
            tf = false;
            if isfield(obj.meta, 'wavelengths') && isnumeric(obj.meta.wavelengths)
                tf = true;
            end
        end
        
        function obj = buffer_bdi(obj, max_free_mem_percentage)
            % attempt to buffer as many ABRDF chunks from a BDI to memory;
            % optionally a fraction can be specified of the available memory
            % that is to be used
            if ~strcmpi(obj.format_str, 'bdi')
                error('buffer_bdi only works with BDIs. non-BDIs are always fully loaded to memory.');
            end
            
            % free memory to ensure the maximum number of available bytes
            % is used
            if obj.is_buffered()
                obj.clear_buffer();
            end
            
            if ~exist('max_free_mem_percentage', 'var')
                max_free_mem_percentage = 0.9;
            end
            
            mem_available = floor(utils.freemem() * max_free_mem_percentage);
            
            mem_required = obj.meta.chunk_size * obj.meta.num_chunks * utils.sizeof('uint16');
            if mem_required > mem_available
                obj.data.num_chunks_in_buffer = floor(obj.meta.num_chunks * mem_available / mem_required);
            else
                obj.data.num_chunks_in_buffer = obj.meta.num_chunks;
            end
            
            if mem_available < obj.meta.chunk_size * utils.sizeof('uint16')
                obj.data.num_chunks_in_buffer = 0;
                error('too little memory available for buffering, one chunk requires %3.2f MB', obj.meta.chunk_size * utils.sizeof('uint16') / 1024 / 1024);
            end
            
            % allocate memory
            obj.data.chunks = zeros(obj.meta.chunk_size, ...
                obj.data.num_chunks_in_buffer, 'uint16');
            
            % read chunks from file
            obj.data.fid = fopen(obj.meta.file_name, 'r');
            for c = 1 : obj.data.num_chunks_in_buffer - 1
                obj.progress(c / obj.data.num_chunks_in_buffer, ...
                    sprintf('(reading chunk %d / %d, %d total)', ...
                    c, obj.data.num_chunks_in_buffer, obj.meta.num_chunks));
                obj.get_bdi_chunk(c);
            end
            c = obj.data.num_chunks_in_buffer;
            % handle last chunk separately because it is usually smaller, but
            % indexing the memory chunks in each of the above iterations would
            % strongly reduce performance
            obj.progress(c / obj.data.num_chunks_in_buffer, ...
                sprintf('(reading chunk %d / %d, %d total)', ...
                c, obj.data.num_chunks_in_buffer, obj.meta.num_chunks));
            obj.get_bdi_chunk(c);
            fclose(obj.data.fid);
            obj.data.chunks_buffered = false(1, obj.meta.num_chunks);
            obj.data.chunks_buffered(1 : obj.data.num_chunks_in_buffer) = true;
            obj.progress();
        end
        
        function obj = clear_buffer(obj)
            % clears BDI buffer to free memory
            if obj.is_bdi()
                obj.data.num_chunks_in_buffer = 0;
                obj.data.chunks_buffered(:) = false;
                obj.data.chunks = [];
            end
        end
        
        function [buffered, num_chunks_buffered, num_chunks_tot] = is_buffered(obj)
            % checks if BDI is buffered
            % (BTFs are always fully buffered, hence the default true)
            buffered = true;
            num_chunks_buffered = 1;
            num_chunks_tot = 1;
            if obj.is_bdi()
                num_chunks_tot = numel(obj.data.chunks_buffered);
                num_chunks_buffered = nnz(obj.data.chunks_buffered);
                buffered = num_chunks_buffered == num_chunks_tot;
            end
        end
        
        function obj = only_use_buffered(obj, value)
            % only use buffered data for a BDI, don't read from file
            % if a bdi is partially buffered, setting this to true will tell the
            % decoder to ignore those pixels that aren't buffered, otherwise, it
            % will (slowly) read those missing pixel from file
            if obj.is_bdi()
                obj.data.only_use_buffered = logical(value(1));
            end
        end
        
        function obj = textures_from_file(obj, value)
            % enable or disable loading textures from file for BDIs
            if obj.is_bdi()
                obj.data.textures_from_file = logical(value(1));
            end
        end
        
        function n = nC(obj)
            % return number of channels
            n = obj.meta.num_channels;
        end
        
        function n = num_channels(obj)
            % return number of channels
            n = obj.meta.num_channels;
        end
        
        function n = nL(obj)
            % return number of unique light directions
            n = obj.meta.nL;
        end
        
        function n = nV(obj)
            % return number of unique view directions
            n = obj.meta.nV;
        end
        
        function n = nLV(obj)
            % return number of unique view directions
            n = obj.meta.nV * obj.meta.nL;
        end
        
        function n = width(obj)
            % return texture width
            n = obj.meta.width;
        end
        
        function n = height(obj)
            % return texture height
            n = obj.meta.height;
        end
        
        function type = data_type(obj)
            % return pixel data type
            type = obj.data.class;
        end
        
        function [L, V] = get_LV_full(obj, linds, vinds)
            % get cartesian product of the sampled light and view directions
            if ~exist('linds', 'var') || isempty(linds) || strcmp(linds, ':')
                linds = 1 : obj.nL;
            end
            
            if ~exist('vinds', 'var') || isempty(vinds) || strcmp(vinds, ':')
                vinds = 1 : obj.nV;
            end
            
            L = obj.get_L_full(linds, vinds);
            V = obj.get_V_full(linds, vinds);
        end
        
        function L = get_L_full(obj, linds, vinds)
            % get light direction for each sample in the data
            if ~exist('linds', 'var') || isempty(linds) || strcmp(linds, ':')
                linds = 1 : obj.nL;
            end
            
            if ~exist('vinds', 'var') || isempty(vinds) || strcmp(vinds, ':')
                vinds = 1 : obj.nV;
            end
            
            L = repmat(obj.meta.L(linds, :), numel(vinds), 1);
        end
        
        function V = get_V_full(obj, linds, vinds)
            % get view direction for each sample in the data
            if ~exist('linds', 'var') || isempty(linds) || strcmp(linds, ':')
                linds = 1 : obj.nL;
            end
            
            if ~exist('vinds', 'var') || isempty(vinds) || strcmp(vinds, ':')
                vinds = 1 : obj.nV;
            end
            
            V = obj.meta.V(vinds, :);
            V = V(:) * ones(1, numel(linds));
            V = reshape(V', [], 3);
        end
        
        function [lin, laz, vin, vaz] = inds_to_angles(obj, l, v)
            % given light and view indices, return the light & view inclination and azimuth angles
            l_sph = utils.cart2sph2(obj.meta.L(l, :));
            v_sph = utils.cart2sph2(obj.meta.V(v, :));
            lin = l_sph(1);
            laz = l_sph(2);
            vin = v_sph(1);
            vaz = v_sph(2);
        end
        
        function [normal_map, height_map] = get_normals_from_height_map(obj)
            % for a BTF resampled to its height field, this function returns the
            % normals computed from this height map
            if ~isfield(obj.data, 'height_map')
                error('This BTF does not store a height field, cannot compute the normals!');
            end
            height_map = obj.data.height_map;
            [height, width] = size(height_map);
            
            ppmm = obj.meta.ppmm;
            dx = 1 / ppmm;
            dy = 1 / ppmm;
            
            [hx, hy] = gradient(height_map);
            normal_map = cat(3, -hx * dy, -hy * dx, ones(height, width) * dx * dy);
            normal_map = normal_map ./ repmat(sqrt(sum(normal_map .^ 2, 3)), [1, 1, 3]);
        end
        
        function name = get_name(obj)
            name = ['BDI ', obj.meta.file_name];
        end
        
        function abrdf = decode_abrdf(obj, x, y, form_cart_prod)
            % decode a full ABRDF for a given pair of x- and y-coordinates
            %
            % no interpolation yet, so x and y should be integers
            % if form_cart_prod is set to true, the cartesian product of
            % all the coordinates in x and y is formed.
            
            if ~exist('form_cart_prod', 'var')
                form_cart_prod = false;
            end
            
            x = uint32(x);
            y = uint32(y);
            
            if ndims(x) ~= ndims(y)
                error('x and y arrays should have the same number of dimensions');
            end
            dim_x = size(x);
            dim_y = size(y);
            
            if ~all(dim_x == dim_y) || form_cart_prod
                % apparently we're supposed to form the cartesian product
                [x, y] = ndgrid(x(:), y(:));
            end
            
            switch lower(obj.format_str)
                case 'bdi'
                    abrdf = obj.decode_bdi_abrdf(x, y);
                case 'dfmf'
                    abrdf = obj.decode_dfmf_abrdf(x, y);
                case 'fmf'
                    abrdf = obj.decode_fmf_abrdf(x, y);
                case 'pvf'
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
                abrdf = bsxfun(@rdivide, abrdf, obj.meta.L(:, 3));
            end
            
            % rearrange as image
            abrdf = reshape(permute(abrdf, [1, 2, 4, 3]), obj.meta.nL, obj.meta.nV, obj.meta.num_channels, []);
        end
        
        function texel = decode_texel(obj, x, y, L, V, form_cart_prod)
            % decode a single texel at given pixel coordinates for given light & view directions
            % 
            % x, y need to be integer indices into the texture space, L and
            % V can either be integers indices in the angular domain, two
            % vectors (either 3D cartesian or 2D inclination-azimuth pairs) or
            % two 3 x N or 2 x N matrices of such vectors.
            % In the latter case, angular interpolation is performed
            % between the color values corresponding to the closest sampled
            % direction pairs.
            % If form_cart_prod is set to true, the cartesian product of
            % all the coordinates in x and y is formed, as well as the cartesian
            % product of all indices in L and V.
            
            if ~exist('form_cart_prod', 'var')
                form_cart_prod = false;
            end
            
            if ndims(x) ~= ndims(y)
                error('x and y arrays should have the same number of dimensions');
            end
            dim_x = size(x);
            dim_y = size(y);
            
            if ndims(L) ~= ndims(V)
                error('x and y arrays should have the same number of dimensions');
            end
            if size(L, 1) ~= 3
                L = L(:, :)';
            end
            if size(V, 1) ~= 3
                V = V(:, :)';
            end
            dim_L = size(L);
            dim_V = size(V);
            
            if ~all(dim_x == dim_y) || form_cart_prod
                % apparently we're supposed to form the cartesian product
                [x, y] = ndgrid(x(:), y(:));
            end
            
            if ~all(dim_L == dim_V) || form_cart_prod
                % apparently we're supposed to form the cartesian product
                L_in = L;
                V_in = V;
                L = zeros(3, dim_L(2), dim_V(2));
                V = zeros(3, dim_L(2), dim_V(2));
                if any(size(L_in) == 1)
                    [L, V] = ndgrid(L_in, V_in);
                elseif any(size(L_in) == 3)
                    for ci = 1 : 3
                        [L(ci, :, :), V(ci, :, :)] = ndgrid(reshape(L_in(ci, :), [], 1), ...
                            reshape(V_in(ci, :), [], 1));
                    end
                end
            end
            
            % L / V could contain indices or actual directions, we need some
            % stupid way to distinguish the two cases...
            if numel(L) == 1 || numel(V) == 1 || ...
                    isinteger(x) || isinteger(y) || ...
                    isinteger(L) || isinteger(V)
                % direct access via indices
                switch lower(obj.format_str)
                    case 'bdi'
                        texel = obj.decode_bdi_texel(x, y, L, V);
                    case 'dfmf'
                        texel = obj.decode_dfmf_texel(x, y, L, V);
                    case 'fmf'
                        texel = obj.decode_fmf_texel(x, y, L, V);
                    case 'pvf'
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
                    texel = bsxfun(@rdivide, texel, obj.meta.L(L, 3));
                end
            elseif any(size(L) == 2) && any(size(V) == 2) || ...
                    any(size(L) == 3) && any(size(V) == 3)
                % access via light and view directions -> interpolate samples
                
                nlv = size(L(:, :), 2);
                assert(nlv == size(V(:, :), 2));
                nxy = numel(x);
                assert(nxy == numel(y));
                
                texel = zeros(nlv, nxy, obj.meta.num_channels, obj.data.class);
                
                % get interpolation indices and weights
                [ls, vs, weights_l, weights_v] = obj.lookup_dirs(L, V);
                for li = 1 : 3
                    l = ls(:, li);
                    for vi = 1 : 3
                        v = vs(:, vi);
                        
                        texel = texel + bsxfun(@times, weights_l(:, li) .* weights_v(:, vi), ...
                            obj.decode_texel(int32(x), int32(y), int32(l), int32(v)));
                    end
                end
            else
                error('error while decoding: either provide readily computed indices, 3D cartesian directions, or two polar angles respectively for light and view.');
            end
        end
        
        function img = decode_texture(obj, L, V, form_cart_prod)
            % decode full texture for given light & view directions
            %
            % L and V can either be integer indices in the angular domain,
            % two vectors (either 3D cartesian or 2D inclination-azimuth pairs)
            % or two 3 x N or 2 x N matrices of such vectors. 
            % In the latter two cases, angular interpolation is performed
            % between the color values corresponding to the closest sampled
            % direction pairs.
            % If form_cart_prod is set to true, the cartesian product of
            % all the coordinates in L and V is formed.
            
            if ~exist('form_cart_prod', 'var')
                form_cart_prod = false;
            end
            
            if ndims(L) ~= ndims(V)
                error('x and y arrays should have the same number of dimensions');
            end
            if size(L, 1) ~= 3
                L = L(:, :)';
            end
            if size(V, 1) ~= 3
                V = V(:, :)';
            end
            dim_L = size(L);
            dim_V = size(V);
            
            if ~all(dim_L == dim_V) || form_cart_prod
                % apparently we're supposed to form the cartesian product
                L_in = L;
                V_in = V;
                L = zeros(3, dim_L(2), dim_V(2));
                V = zeros(3, dim_L(2), dim_V(2));
                for ci = 1 : 3
                    [L(ci, :, :), V(ci, :, :)] = ndgrid(reshape(L_in(ci, :), [], 1), ...
                        reshape(V_in(ci, :), [], 1));
                end
            end
            
            % do L / V contain integer indices or vectors?
            if numel(L) ==  1 || numel(V) == 1 || isinteger(L) || isinteger(V)
                L = round(L);
                V = round(V);
                % direct access via indices
                switch lower(obj.format_str)
                    case 'bdi'
                        img = obj.decode_bdi_texture(L, V);
                    case 'dfmf'
                        img = obj.decode_dfmf_texture(L, V);
                    case 'fmf'
                        img = obj.decode_fmf_texture(L, V);
                    case 'pvf'
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
                    img = bsxfun(@rdivide, img, reshape(obj.meta.L(L, 3), 1, 1, []));
                end
                
                % rearrange as image
                img = permute(reshape(img, obj.meta.height, obj.meta.width, [], ...
                    obj.meta.num_channels), [1, 2, 4, 3]);
            elseif any(size(L) == 2) || any(size(L) == 3)
                % access via light and view directions -> interpolate samples
                nlv = size(L(:, :), 2);
                assert(nlv == size(V(:, :), 2));
                
                img = zeros(obj.meta.height, obj.meta.width, ...
                    obj.meta.num_channels, nlv, obj.data.class);
                [ls, vs, weights_l, weights_v] = obj.lookup_dirs(L, V);
                weights_l = reshape(weights_l', 1, 1, 3, []);
                weights_v = reshape(weights_v', 1, 1, 3, []);
                for li = 1 : 3
                    l = ls(:, li);
                    for vi = 1 : 3
                        v = vs(:, vi);
                        
                        img = img + bsxfun(@times, weights_l(:, :, li, :) .* weights_v(:, :, vi, :), ...
                            obj.decode_texture(int32(l), int32(v)));
                    end
                end
            else
                error('error while decoding: either provide readily computed indices, 3D cartesian directions, or two polar angles respectively for light and view.');
            end
        end
        
        function textures = get_eigen_textures(obj)
            % return the eigen textures of a matrix factorization-based
            % compressed BTF
            switch lower(obj.format_str)
                case 'dfmf'
                    textures = obj.get_dfmf_eigen_textures();
                case 'fmf'
                    textures = obj.get_fmf_eigen_textures();
                case 'pvf'
                    textures = obj.get_pvf_eigen_textures();
                otherwise
                    error('extraction of eigen textures is only possible with (D)FMF BTFs or PVF BTFs.');
            end
        end
        
        function abrdfs = get_eigen_abrdfs(obj)
            % return the eigen ABRDFs of a matrix factorization-based
            % compressed BTF
            switch lower(obj.format_str)
                case 'dfmf'
                    abrdfs = obj.get_dfmf_eigen_abrdfs();
                case 'fmf'
                    abrdfs = obj.get_fmf_eigen_abrdfs();
                case 'pvf'
                    abrdfs = obj.get_pvf_eigen_abrdfs();
                otherwise
                    error('extraction of eigen ABRDFs is only possible with (D)FMF BTFs or PVF BTFs.');
            end
        end
        
        function tensor = get_5d_tensor(obj)
            % extract the full data tensor as a 5-dimensional array (num_lights
            % x num_views x width x height x num_channels)
            switch lower(obj.format_str)
                case 'bdi'
                    if numel(obj.data.chunks) == obj.meta.num_channels * ...
                            obj.meta.nL * obj.meta.nV * obj.meta.width * obj.meta.height
                        tensor = halfprecision(reshape(obj.data.chunks, obj.meta.num_channels, ...
                            obj.meta.nL, obj.meta.nV, obj.meta.width, obj.meta.height), 'single');
                    else
                        error(['Tensor can only be returned for fully buffered BDI. ', ...
                            'Please call buffer_bdi().']);
                    end
                case 'dfmf'
                    tensor = cell(obj.meta.num_channels, 1);
                    for c = 1 : obj.meta.num_channels
                        tensor{c} = obj.data.U{c} * obj.data.SxV{c}';
                    end
                    tensor = cat(3, tensor{:});
                    tensor = obj.undecorrelate(tensor);
                    tensor = reshape(tensor, obj.meta.nL, obj.meta.nV, ...
                        obj.meta.width, obj.meta.height, obj.meta.num_channels);
                case 'fmf'
                    error('not implemented yet!');
                case 'pvf'
                    error('not implemented yet!');
            end
        end
        
        function obj = set_verbose(obj, verbose)
            % toggle object's verbosity;
            % at the moment this only enables or disables calls to the object's
            % status callback
            obj.verbose = logical(verbose(1));
        end
        
        function obj = set_progress_callback(obj, progress_callback)
            % set a progress callback function, that is called during lengthy processes;
            % the handle needs to be to a function that takes two parameters,
            % the first one is a scalar between 0 and 1, indicating a
            % percentage, the second one is an optional string containing a
            % short description of the process; the function also needs to be
            % able to be called without any arguments, which can be used to hide
            % or reset gui elements related to the progress display
            obj.progress_fcn = progress_callback;
        end
        
        function obj = set_cosine_flag(obj, flag)
            % set the boolean flag that indicates whether the cosine of
            % theta_light should be divided out of the data or not
            obj.meta.cosine_flag = logical(flag(1));
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
            L = utils.add_bottom_ring(L');
            V = utils.add_bottom_ring(V');
            
            % set up Delaunay trianglulation of parabolic map of both light
            % and view directions
            L_par = utils.cart2par(L);
            V_par = utils.cart2par(V);
            obj.TriL = delaunayTriangulation(L_par');
            obj.TriV = delaunayTriangulation(V_par');
        end
        
        function [L_idxs, V_idxs, L_weights, V_weights] = lookup_dirs(obj, L, V)
            % determine indices and weights for angular interpolation
            
            % transform input directions to parabolic coordinates
            if size(L, 1) == 3 && size(V, 1) == 3
                L_par = utils.cart2par(L);
                V_par = utils.cart2par(V);
            elseif size(L, 1) == 2 && size(V, 1) == 2
                % we assume L and V to be in parabolic coordinates already!
                % unfortunately, there is no sane way to distinguish those
                % from spherical coordinates
                L_par = L;
                V_par = V;
            end
            
            % lookup in Delaunay triangulation for both directions
            [L_ti, L_weights] = obj.TriL.pointLocation(reshape(L_par, 2, [])');
            [V_ti, V_weights] = obj.TriV.pointLocation(reshape(V_par, 2, [])');
            valid = ~isnan(L_ti) & ~isnan(V_ti);
            
            % we cannot return indices into the bottom ring, as we don't
            % have any data available there
            L_idxs = (obj.meta.nL + 1) * ones(numel(valid), 3);
            V_idxs = (obj.meta.nV + 1) * ones(numel(valid), 3);
            L_idxs(valid, :) = obj.TriL.ConnectivityList(L_ti(valid), :);
            V_idxs(valid, :) = obj.TriV.ConnectivityList(V_ti(valid), :);
            
            % hence, we're only interested in those indices and weights
            % that are actually containted in the sampling
            L_weights(L_idxs > obj.meta.nL) = 0;
            V_weights(V_idxs > obj.meta.nV) = 0;
            L_idxs(L_idxs > obj.meta.nL) = 1;
            V_idxs(V_idxs > obj.meta.nV) = 1;
        end
        
        function progress(obj, value, str, varargin)
            % this function is called by members of btf objects to send status updates
            if obj.verbose
                try %#ok<TRYNC>
                    if exist('value', 'var')
                        obj.progress_fcn(value, str, varargin{:});
                    else
                        obj.progress_fcn();
                    end
                end
            end
        end
        
        function default_progress_fcn(obj, value, str) %#ok<INUSL>
            % this is just a very simple function for displaying progress
            % updates on the matlab prompt
            if ~exist('str', 'var')
                str = '';
            end
            
            if exist('value', 'var')
                fprintf('%03.2f%% (%s)...\n', 100 * value, str);
            end
        end
    end
end
