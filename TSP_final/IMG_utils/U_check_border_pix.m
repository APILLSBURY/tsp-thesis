function border = U_check_border_pix(IMG_label, index, cur_label)
    [xdim, ydim] = size(IMG_label);
    [x, y] = get_x_and_y_from_index(index, xdim);
    x1 = max(1, x-1);
    y1 = max(1, y-1);
    x2 = min(xdim, x+1);
    y2 = min(ydim, y+1);
    
    %check to see if a label was specified
    if nargin < 3
        cur_label = IMG_label(x, y);
    end

    %check to see if any of the adjacent pixels are different than
    %the current one
    border = (cur_label ~= IMG_label(x1, y) || cur_label ~= IMG_label(x, y1) || ...
        cur_label ~= IMG_label(x2, y) || cur_label ~= IMG_label(x, y2));
end