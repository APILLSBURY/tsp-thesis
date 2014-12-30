% --------------------------------------------------------------------------
% -- fix_borders
% --   Fixes the border linked list and the border_ptr image for a single
% -- super pixel.
% --
% --   parameters:
% --     - label : the label image
% --     - border_ptr : the border_ptr image
% --     - (xdim, ydim) : the size of the image
% --------------------------------------------------------------------------
function SP = SP_fix_borders(SP, label)
    [xdim, ydim] = size(label);
    % update the border pixels
    SP.borders = false(size(SP.borders));
    for index=find(SP.pixels)'
        [x, y] = get_x_and_y_from_index(index, xdim);
        cur_label = label(x, y);
        SP.borders(index) = (x>1 && label(x-1, y)~=cur_label) || (y>1 && label(x, y-1)~=cur_label) || ...
                            (x<xdim && label(x+1, y)~=cur_label) || (y<ydim && label(x, y+1)~=cur_label);
    end
end