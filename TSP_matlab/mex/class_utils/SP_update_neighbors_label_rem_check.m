% --------------------------------------------------------------------------
% -- update_neighbors_label_rem_check
% --   Checks to see if the neighbor count at index should be decremented
% -- by removing one neighbor of label neighbor_label. If so, it decrements.
% -- The neighboring label should be changed before calling this function.
% --
% --   parameters:
% --     - label : the label image
% --     - index : the index bordering the removed pixel
% --     - (xdim,ydim) : dimensions of image
% --     - neighbor_label : the label of the neighbor to decrement
% --------------------------------------------------------------------------
function SP = SP_update_neighbors_label_rem_check(SP, label, x, y, neighbor_label)
    % check to see how many of the neighbors have label neighbor_label
    [xdim, ydim] = size(label);

    if ~((x>1 && label(x-1, y)==neighbor_label) || (y>1 && label(x, y-1)==neighbor_label) || ...
            (x<xdim && label(x+1, y)==neighbor_label) || (y<ydim && label(x, y+1)==neighbor_label))
        SP = SP_update_neighbors_label_rem(SP, neighbor_label);
    end
end