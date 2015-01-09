% --------------------------------------------------------------------------
% -- U_update_border_changed
% --   Updates all border linked lists and pointers for the neighbors of a
% -- changed pixel at index.
% --
% --   parameters:
% --     - index : the index of the recently changed pixel
% --------------------------------------------------------------------------
function IMG_SP = U_update_border_changed(IMG_label, IMG_SP, index)
    [xdim, ydim] = size(IMG_label);
    [x, y] = get_x_and_y_from_index(index, xdim);
    if (x>1)
        curLabel = IMG_label(x-1, y);
        if curLabel > 0
            IMG_SP(curLabel).borders(index-1) = U_check_border_pix(IMG_label, index-1);
        end
    end
    if (y>1)
        curLabel = IMG_label(x, y-1);
        if curLabel > 0
            IMG_SP(curLabel).borders(index-xdim) = U_check_border_pix(IMG_label, index-xdim);
        end
    end
    if (x<xdim)
        curLabel = IMG_label(x+1, y);
        if curLabel > 0
            IMG_SP(curLabel).borders(index+1) = U_check_border_pix(IMG_label, index+1);
        end
    end
    if (y<ydim)
        curLabel = IMG_label(x, y+1);
        if curLabel > 0
            IMG_SP(curLabel).borders(index+xdim) = U_check_border_pix(IMG_label, index+xdim);
        end
    end
end