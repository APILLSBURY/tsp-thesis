% --------------------------------------------------------------------------
% -- update_neighbors_add_self
% --   Updates the neighbor lists and counts by adding one particular pixel
% -- at index. Does not update the neighboring neighbor lists.
% --
% --   parameters:
% --     - label : the label image
% --     - index : the index of the added pixel
% --     - (xdim,ydim) : dimensions of image
% --------------------------------------------------------------------------
function SP = SP_update_neighbors_add_self(SP, label, index)
    [xdim, ydim] = size(label);
    [x, y] = get_x_and_y_from_index(index, xdim);
    cur_label = label(x, y);

    % don't do all neighbors, only the unique ones
    if x>1
        llabel = label(x-1, y);
    else
        llabel = 0;
    end
    if y>1
        ulabel = label(x, y-1);
    else
        ulabel = 0;
    end
    if x<xdim
        rlabel = label(x+1, y);
    else
        rlabel = 0;
    end
    if y<ydim
        dlabel = label(x, y+1);
    else
        dlabel = 0;
    end

    if (llabel>0 && llabel~=cur_label)
        SP = SP_update_neighbors_label_add(SP, llabel);
    end
    if (ulabel>0 && ulabel~=cur_label)
        SP = SP_update_neighbors_label_add(SP, ulabel);
    end
    if (rlabel>0 && rlabel~=cur_label)
        SP = SP_update_neighbors_label_add(SP ,rlabel);
    end
    if (dlabel>0 && dlabel~=cur_label)
        SP = SP_update_neighbors_label_add(SP, dlabel);
    end
end