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
    neighbor_labels = zeros(1, 4);
    if x>1
        neighbor_labels(1) = label(x-1, y);
    end
    if y>1
        neighbor_labels(2) = label(x, y-1);
    end
    if x<xdim
        neighbor_labels(3) = label(x+1, y);
    end
    if y<ydim
        neighbor_labels(4) = label(x, y+1);
    end
    
    for neighbor_label=unique(neighbor_labels)
        if neighbor_label>0 && neighbor_label~=cur_label
            SP.neighbors(neighbor_label) = SP.neighbors(neighbor_label)+1;
        end
    end
end