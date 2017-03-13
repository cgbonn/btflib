% *************************************************************************
% * Copyright 2017 University of Bonn
% *
% * authors:
% *  - Sebastian Merzbach <merzbach@cs.uni-bonn.de>
% *
% * file creation date: 2017-01-25
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
% simple unit tests for checking the functionality of the OpenEXR reading and
% writing functions. run by calling runtests('exr.test');
rng(0);

% set up single precision floating point image
im = rand(300, 200, 33, 'single');
wls = utils.sprintf2('%3.2f', linspace(400, 700, size(im, 3)));

% set up unsigned integer image
imui = randi(255, 300, 200, size(im, 3), 'uint32');

%% write and read a simple RGB floating point precision image
imrgb = rand(10, 10, 3, 'single');
exr.write(imrgb, fullfile(tempdir, 'test_rgb.exr'), 'single');
[imrrgb_outputfile, wlsrrgb] = exr.read(fullfile(tempdir, 'test_rgb.exr'), true);
assert(all(imrrgb_outputfile(:) == imrgb(:)), ...
    'difference to single precision RGB file is nonzero!');

%% write to single precision from single precision image
exr.write(im, fullfile(tempdir, 'test_single.exr'), 'single', wls);
[imrs, wlsr] = exr.read(fullfile(tempdir, 'test_single.exr'));
assert(all(im(:) == imrs(:)), ...
    'difference between single precision image after writing and reading is nonzero!');

%% write to half precision from single precision array
exr.write(im, fullfile(tempdir, 'test_half_from_single.exr'), 'half', wls);
imrhs = exr.read(fullfile(tempdir, 'test_half_from_single.exr'));
assert(mean(abs(im(:) - halfprecision(imrhs(:), 'single'))) < 3e-4, ...
    'difference between half precision image after writing and reading is too big!');
% that's what it's supposed to be like, but apparently OpenEXR's half.h behaves
% differently than the one used in halfprecision.c, which leads to slightly
% different results when converting from floats...
% assert(all(halfprecision(im(:)) == imrhs(:)), ...
%     'difference between half precision image after writing and reading is nonzero!');

%% write to half precision from half precision (uint16) array
exr.write(halfprecision(im), fullfile(tempdir, 'test_half_from_half.exr'), 'half', wls);
imrhh = exr.read(fullfile(tempdir, 'test_half_from_half.exr'));
assert(all(halfprecision(im(:)) == imrhh(:)), ...
    'difference between half precision image after writing and reading is nonzero!');

%% write to UINT file from single precision array
exr.write(im, fullfile(tempdir, 'test_uint_from_single.exr'), 'uint', wls);
imruis = exr.read(fullfile(tempdir, 'test_uint_from_single.exr'), false);
assert(all(uint32(im(:) - 0.5) == imruis(:)), ...
    'difference between unsigned int image after writing and reading is nonzero!');

%% write to UINT file from unsigned integer array
exr.write(imui, fullfile(tempdir, 'test_uint.exr'), 'uint', wls);
imrui = exr.read(fullfile(tempdir, 'test_uint.exr'), false);
assert(all(imui(:) == imrui(:)), ...
    'difference between unsigned int image after writing and reading is nonzero!');
