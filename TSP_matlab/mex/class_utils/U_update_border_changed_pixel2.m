% --------------------------------------------------------------------------
% -- U_update_border_changed_pixel
% --   Updates all border linked lists and pointers at index
% --
% --   parameters:
% --     - index : the index of the recently changed pixel
% --------------------------------------------------------------------------
function IMG = U_update_border_changed_pixel2(IMG, index)
    [x, y] = get_x_and_y_from_index(index, IMG.xdim);
    curLabel = IMG.label(x, y);
    if curLabel>0
        IMG.SP(curLabel).borders(index) = U_check_border_pix(IMG, index);
    end
end