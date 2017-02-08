% *************************************************************************
% * Copyright 2017 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * last modification date: 2017-02-03
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
% This function reads in a folder with rectified HDR images to create a
% BDI. Currently, the BDI is fully buffered in memory, so beware memory
% limits. It is possible to read images selectively by specifying a regular
% expression that filters out certain combinations of light and view
% angles. Furthermore, to save memory one can specify a region of interest.
%
% Currently, only images in the format of the SPECTRAL dataset from
% University of Bonn's BTFDBB are supported. More formats will be added in
% the future.
%
% Usage: btf_obj = imgs2bdi(folder, filter, roi, output_format, verbose),
% where
% - folder is the path to the folder containing the images
% - filter is a string with a regular expression for selecting only
%   specific images; if empty, it defaults to 'tv.*', which selects all
%   images; an example for selecting only specific angles:
%   ''tv..5.0_pv..0.0_tl..5.0_pl..0.0.*' selects the theta angles to end
%   in 5.0, and the phi angles to end in 0.0; one example for the file
%   names is: tv045.0_pv020.0_tl015.0_pl120.0.exr
% - roi is a 2 x 3 array which is interpreted as follows: roi = [x_min,
%   x_stride, x_max; y_min, y_stride; y_max]; it defaults to the full
%   texture resolution with a stride of 1 in both x and y
% - output_format allows to select the data type the buffered images are
%   stored in; possible values are 'half' (default) or 'single'; half
%   stores the data with 2 byte half precision floats (internally as
%   uint16, which is converted to single precision after evaluating the
%   BTF)
% - verbose is an integer from 0 to 3, which allows to enable status
%   information
%
% This function makes 
function btf_obj = imgs2bdi(folder, filter, roi, output_format, verbose)
    if ~exist('filter', 'var') || isempty(filter)
        filter = 'tv.*';
    end
    
    if ~exist('roi', 'var') || isempty(roi)
        roi = [];
    end
    
    if ~exist('output_format', 'var')
        output_format = 'half';
    end
    
    assert(ismember(output_format, {'half', 'single'}), ...
        'output_format must be ''half'' or ''single''.');
    
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 2;
    end
    
    folder = tb.make_absolute(folder);
    
    % get file names and sort them by light / view angles
    fnames = dir(folder);
    folders = repmat({folder}, numel(fnames), 1);
    fnames = {fnames.name}';
    matching = regexpi(fnames, filter);
    matching = cellfun(@(x) ~isempty(x), matching);
    folders = folders(matching);
    fnames = fnames(matching);
    fnames = cellfun(@(x, y) fullfile(x, y), folders, fnames, 'UniformOutput', false);
    
    % check extension
    [~, ~, exts] = cellfun(@(x) fileparts(x), fnames, 'UniformOutput', false);
    assert(all(cellfun(@(x) strcmpi(x, exts{1}), exts)), 'images are not all in the same format!');
    ext = exts{1}(2 : end);
    
    % sort by angles
    tokens = regexpi(fnames, 'tv(\d+\.?\d*)_pv(\d+\.?\d*)_tl(\d+\.?\d*)_pl(\d+\.?\d*)', 'tokens');
    [tv, pv, tl, pl] = cellfun(@(x) deal(str2double(x{1}{1}), str2double(x{1}{2}), ...
        str2double(x{1}{3}), str2double(x{1}{4})), tokens);
    angles = [tv, pv, tl, pl];
    [angles, perm] = sortrows(angles);
    angles = angles ./ 180 * pi;
    [tv, pv, tl, pl] = deal(angles(:, 1), angles(:, 2), angles(:, 3), angles(:, 4));
    fnames = fnames(perm);
    
    % compute light and view in cartesian coordinates
    nlv = numel(fnames);
    tpl = unique([tl, pl], 'rows');
    tpv = unique([tv, pv], 'rows');
    L = utils.sph2cart2(tpl);
    V = utils.sph2cart2(tpv);
    nl = size(tpl, 1);
    nv = size(tpv, 1);
    
    assert(nl * nv == nlv);
    
    if ~strcmpi(ext, 'exr')
        error('currently image-based BTFs are only supported in OpenEXR format.');
    end
    
    meta = exr.query(fnames{1});
    h = meta.height;
    w = meta.width;
    nc = meta.num_channels;
    num_scan_lines_per_chunk = 4;
    
    % deal with region of interest
    if isempty(roi)
        roi = [1, 1, w; 1, 1, h];
    end
    if roi(1, 3) < 1
        roi(1, 3) = w;
    end
    if roi(2, 3) < 1
        roi(2, 3) = h;
    end
    
    if any(roi(:, 3) > [w; h]) || any(roi(:, 1) < 1)
        warning('Region of interest exceeds texture dimensions, will be clamped!');
    end
    roi(1, 1) = max(1, roi(1, 1));
    roi(1, 3) = min(w, roi(1, 3));
    roi(2, 1) = max(1, roi(2, 1));
    roi(2, 3) = min(h, roi(2, 3));
    
    % we don't apply the stride yet because we will use imresize instead to
    % get nicer downsampling results
    stride_x = roi(1, 2);
    stride_y = roi(2, 2);
    xs = roi(1, 1) : roi(1, 3);
    ys = roi(2, 1) : roi(2, 3);
    
    w = numel(xs(1 : stride_x : end));
    h = numel(ys(1 : stride_y : end));
    
    % pre-allocate memory
    wlsstrs = cell(1, nlv);
    wls = cell(1, nlv);
    chunk_format = 'single';
    if strcmpi(output_format, 'half')
        chunk_format = 'uint16';
    end
    chunks = zeros(nc, nl * nv, w, h, chunk_format);
    
    if verbose > 1 && exist('multiWaitbar', 'file') == 3
        multiWaitbar('reading images');
    end
    for fi = 1 : nlv
        if verbose > 0
            fprintf('reading image %d / %d (%3.2f%%)\n', fi, nlv, 100 * fi / nlv);
            if verbose > 1 && exist('multiWaitbar', 'file') == 3
                multiWaitbar('reading images', fi / nlv);
            end
        end
        
        % read image to single precision float array
        [im, wlsstrs{fi}] = exr.read(fnames{fi}, true);
        
        % resize image (image should be in float format for that purpose)
        im = imresize(im(ys, xs, :), [h, w], 'bilinear');
        
        % convert to half if requested
        if isa(im, 'single') && strcmp(output_format, 'half')
            im = halfprecision(im);
        end
        
        % store image in chunk buffer
        chunks(:, fi, :, :) = permute(im, [3, 2, 1]);
        
        wls{fi} = str2double(wlsstrs{fi});
    end
    if verbose > 1 && exist('multiWaitbar', 'file') == 3
        multiWaitbar('reading images', 'Close');
    end
    
    % check wavelengths are the same for all images
    assert(all(cellfun(@(x) all(x == wls{1}), wls)), ...
        'wavelength sampling is not the same for all images!');
    
    % create BDI from images
    chunks = reshape(chunks, nc * nl * nv, w * h);
    meta = utils.create_meta_struct('num_lights', nl, 'num_views', nv, ...
        'num_channels', nc, 'width', w, 'height', h, ...
        'light_dirs', L, 'view_dirs', V, 'cosine_flag', false, ...
        'wavelengths', wls{1}, 'channel_names', wlsstrs{1}, ...
        'num_scan_lines_per_chunk', num_scan_lines_per_chunk);
    data = struct('chunks', chunks);
    btf_obj = btf('BDI', meta, data);

    if verbose > 2
        btfview(b);
    end
end
