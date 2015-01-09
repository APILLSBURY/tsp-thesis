% --------------------------------------------------------------------------
% -- U_update_neighbors_rem
% --   Updates all neighbor lists when a pixel is "removed". A subsequent
% -- call to update_neighbors_add should be completed right after this one.
% -- The neighboring label should be changed before calling this function.
% --
% --   parameters:
% --     - old_label : the label of the pixel before it was changed
% --     - index : the index bordering the removed pixel
% --------------------------------------------------------------------------
function IMG_SP = U_update_neighbors_rem(IMG_label, IMG_SP, old_label, index)
    [xdim, ydim] = size(IMG_label);
    if old_label>0
        [x, y] = get_x_and_y_from_index(index, xdim);
        neighbor_labels = zeros(1, 4);
        if x>1
            neighbor_labels(1) = IMG_label(x-1, y);
        end
        if x<xdim
            neighbor_labels(2) = IMG_label(x+1, y);
        end
        if y>1
            neighbor_labels(3) = IMG_label(x, y-1);
        end
        if y<ydim
            neighbor_labels(4) = IMG_label(x, y+1);
        end

        unique_neighbors = unique(neighbor_labels);
        for nindex = 1:length(unique_neighbors);
            neighbor_label=unique_neighbors(nindex);
            if neighbor_label > 0 && neighbor_label ~= old_label
                num_neighbors = IMG_SP(old_label).neighbors(neighbor_label);
                if num_neighbors > 0
                    IMG_SP(old_label).neighbors(neighbor_label) = num_neighbors-1;
                else
                    disp('U_update_neighbors_rem Trying to remove a neighbor that was never added');
                end
            end
        end

        % update the neighbors' neighbors list. (not a typo!)
        if neighbor_labels(1)>0 && neighbor_labels(1)~=old_label
            IMG_SP(neighbor_labels(1)) = SP_update_neighbors_label_rem_check(IMG_SP(neighbor_labels(1)), IMG_label, x-1, y, old_label);
        end
        if neighbor_labels(2)>0 && neighbor_labels(2)~=old_label
            IMG_SP(neighbor_labels(2)) = SP_update_neighbors_label_rem_check(IMG_SP(neighbor_labels(2)), IMG_label, x+1, y, old_label);
        end
        if neighbor_labels(3)>0 && neighbor_labels(3)~=old_label
            IMG_SP(neighbor_labels(3)) = SP_update_neighbors_label_rem_check(IMG_SP(neighbor_labels(3)), IMG_label, x, y-1, old_label);
        end
        if neighbor_labels(4)>0 && neighbor_labels(4)~=old_label
            IMG_SP(neighbor_labels(4)) = SP_update_neighbors_label_rem_check(IMG_SP(neighbor_labels(4)), IMG_label, x, y+1, old_label);
        end
    end
end