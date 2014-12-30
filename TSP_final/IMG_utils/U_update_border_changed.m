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
        curLabel = IMG.label(x-1, y);
        if curLabel > 0
            IMG.SP(curLabel).borders(index-1) = U_check_border_pix(IMG, index-1);
        end
    end
    if (y>1)
        curLabel = IMG.label(x, y-1);
        if curLabel > 0
            IMG.SP(curLabel).borders(index-IMG.xdim) = U_check_border_pix(IMG, index-IMG.xdim);
        end
    end
    if (x<IMG.xdim)
        curLabel = IMG.label(x+1, y);
        if curLabel > 0
            IMG.SP(curLabel).borders(index+1) = U_check_border_pix(IMG, index+1);
        end
    end
    if (y<IMG.ydim)
        curLabel = IMG.label(x, y+1);
        if curLabel > 0
            IMG.SP(curLabel).borders(index+IMG.xdim) = U_check_border_pix(IMG, index+IMG.xdim);
        end
    end
end