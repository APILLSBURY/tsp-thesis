% --------------------------------------------------------------------------
% -- U_update_border_changed
% --   Updates all border linked lists and pointers for the neighbors of a
% -- changed pixel at index.
% --
% --   parameters:
% --     - index : the index of the recently changed pixel
% --------------------------------------------------------------------------
function IMG = U_update_border_changed(IMG, index)
    [x, y] = get_x_and_y_from_index(index, IMG.xdim);
    if (x>1)
        IMG = U_update_border_changed_pixel(IMG, index-1);
    end
    if (y>1)
        IMG = U_update_border_changed_pixel(IMG, index-IMG.xdim);
    end
    if (x<IMG.xdim)
        IMG = U_update_border_changed_pixel(IMG, index+1);
    end
    if (y<IMG.ydim)
        IMG = U_update_border_changed_pixel(IMG, index+IMG.xdim);
    end
end