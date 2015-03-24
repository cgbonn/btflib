function obj = crop_dfmf(obj, roi, strides)
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    
    inds_x = roi(1, 1) : strides(1) : roi(1, 2);
    inds_y = roi(2, 1) : strides(2) : roi(2, 2);
    
    [inds_x, inds_y] = meshgrid(inds_x, inds_y);
    linds_xy = sub2ind([obj.meta.width, obj.meta.height], inds_x', inds_y');
    
    % determine transposition of BTF-matrix
    if size(obj.data.U{1}, 1) == nL * nV
        % images stacked row-wise
        for c = 1 : obj.meta.num_channels
            obj.data.SxV{c} = obj.data.SxV{c}(linds_xy, :);
        end
    elseif size(obj.data.U, 1) == h * w
        % images stacked column-wise
        for c = 1 : obj.meta.num_channels
            obj.data.U{c} = obj.data.U{c}(linds_xy, :);
        end
    else
        error('unknown format of BTF U-component');
    end
end