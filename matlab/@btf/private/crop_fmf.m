function obj = crop_fmf(obj, roi, strides)
    nC = obj.meta.num_channels;
    nL = obj.meta.nL;
    nV = obj.meta.nV;
    h = obj.meta.height;
    w = obj.meta.width;
    
    inds_x = roi(1, 1) : strides(1) : roi(1, 2);
    inds_y = roi(2, 1) : strides(2) : roi(2, 2);
    [inds_x, inds_y] = meshgrid(inds_x, inds_y);
    
    % determine transposition of BTF-matrix
    if size(obj.data.U,1) == nC * nL * nV
        % images stacked row-wise
        linds_xy = sub2ind([w, h], inds_x', inds_y');
        obj.data.SxV = obj.data.SxV(linds_xy, :);
    elseif size(obj.data.U,1) == nC * h * w
        % images stacked column-wise
        inds_x = inds_x';
        inds_y = inds_y';
        linds_cxy = sub2ind([nC, w, h], repmat(1 : nC, numel(inds_x), 1), ...
            repmat(inds_x(:), 1, nC), repmat(inds_y(:), 1, nC));
        obj.data.U = obj.data.U(linds_cxy', :);
    else
        error('unknown format of BTF U-component');
    end    
end