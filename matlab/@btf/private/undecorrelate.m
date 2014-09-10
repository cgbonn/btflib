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
% In the decorrelated full matrix factorization (DFMF) format the BTF's
% color values are decorrelated in one of many alternative formats (e.g.
% the YUV color space). Each of the resulting decorrelated color channels
% is then arranged as a huge matrix which is factorized with an SVD. With
% this method, better compression ratios can be obtained by storing fewer
% components for some of the channels.
% This method restores the original RGB color values from the decorrelated
% channels.
function values = undecorrelate(obj, values, varargin)
    if isfield(obj.meta, 'drr_eps')
        drr_eps = obj.meta.drr_eps;
        if isfield(obj.meta, 'drr_offset')
            drr_eps = drr_eps + obj.meta.drr_offset;
        end
    else
        drr_eps = 1e-5;
    end
    
    if isfield(obj.meta, 'color_model')
        method = obj.meta.color_model;
    else
        method = 0;
    end

    mat_YUV2RGB = [ 1, 0, 1.13983; 1, -0.39465, -0.5806; 1, 2.03211, 0];
    mat_YCoCg2RGB = [1, 1, -1; 1, 0, 1; 1, -1, -1];

    % rearrange data to allow application of undecorrelation matrices
    values = reshape(values, [], obj.meta.num_channels);
    switch method
        case 0 % SIMPLE aka RGB
        case 1 % YUV
            values = values * mat_YUV2RGB';
        case 2 % LAB
            C = makecform('lab2srgb');
            values = applycform(values, C);
        case 3 % PCA
            mat_color_trafo = obj.meta.color_transformation_matrix;
            values = values * mat_color_trafo;
            values = values + repmat(obj.meta.color_mean(:)', size(values, 1), 1);
        case 4 % YCoCg
            values = values * mat_YCoCg2RGB';
        case 5 % SPECTRAL
            first = 1;
            f_sum = values(:,1);
            for c = 2 : size(values, 2)
                first = first - values(:, c);
                values(:, c) = values(:, c) .* f_sum;
            end
            values(:, 1) = first .* f_sum;
        case 6 % SPECTRAL_LINEAR
            f_sum = values(:,1);
            div = 1 / size(values, 2);
            f_sum = f_sum ./ div;
            values(:, 1) = f_sum;
            for c = 2 : size(values, 2)
                values(:, 1) = values(:, 1) - values(:, c) * div;
                values(:, c) = f_sum + values(:, c) .* f_sum;
            end
        case 9 % Ylog_UV
            values(:,1) = exp(values(:,1)) - drr_eps;
            values = values * mat_YUV2RGB';
        case 10 % Ylog_Ulog_Vlog
            values(:,1) = exp(values(:,1)) - drr_eps;
            values(:,2:end) = sign(values(:,2:end)) .* exp(abs(values(:,2:end)) + log(drr_eps));
            values = values * mat_YUV2RGB';
        case 11 % Ylog_UdivY_VdivY
            values(:,1) = exp(values(:,1)) - drr_eps;
            values(:,2) = values(:,2) .* (values(:,1) + drr_eps);
            values(:,3) = values(:,3) .* (values(:,1) + drr_eps);
            values = values * mat_YUV2RGB';
        case 12 % log(RGB) aka SIMPLE_LOG
            values = exp(values) - drr_eps;
        otherwise
            error( 'unsupported color model' );
    end
end
