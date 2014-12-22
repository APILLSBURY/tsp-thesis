% --------------------------------------------------------------------------
% -- update_neighbors_label_rem_check
% --   Checks to see if the neighbor count at index should be decremented
% -- by removing one neighbor of label neighbor_label. If so, it decrements.
% -- The neighboring label should be changed before calling this function.
% -- The count is decremented if there are no neighbors with neighbor_label.
% --
% --   parameters:
% --     - label : the label image
% --     - index : the index bordering the removed pixel
% --     - (xdim,ydim) : dimensions of image
% --     - neighbor_label : the label of the neighbor to decrement
% --------------------------------------------------------------------------
function SP = SP_update_neighbors_label_rem_check(SP, label, x, y, neighbor_label)
    % if no neighbors of this pixel have label neighbor_label, subtract it
    % from neighbors
    [xdim, ydim] = size(label);

    if ~((x>1 && label(x-1, y)==neighbor_label) || (y>1 && label(x, y-1)==neighbor_label) || ...
            (x<xdim && label(x+1, y)==neighbor_label) || (y<ydim && label(x, y+1)==neighbor_label))
        neighbor_val = SP.neighbors(neighbor_label);
        if neighbor_val > 0
            SP.neighbors(neighbor_label) = neighbor_val - 1;
        else
            disp('Trying to remove a neighbor that was never added');
        end
    end
end