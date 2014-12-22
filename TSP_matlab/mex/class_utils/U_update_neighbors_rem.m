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
function IMG = U_update_neighbors_rem(IMG, old_label, index)
    if (old_label>0)
        [x, y] = get_x_and_y_from_index(index, IMG.xdim);
        neighbor_labels = ones(1, 4) * old_label;
        if x>1
            neighbor_labels(1) = IMG.label(x-1, y);
        end
        if x<IMG.xdim
            neighbor_labels(2) = IMG.label(x+1, y);
        end
        if y>1
            neighbor_labels(3) = IMG.label(x, y-1);
        end
        if y<IMG.ydim
            neighbor_labels(4) = IMG.label(x, y+1);
        end

        for neighbor_label=unique(neighbor_labels)
            if neighbor_label > 0 && neighbor_label ~= old_label
                num_neighbors = IMG.SP(old_label).neighbors(neighbor_label);
                if num_neighbors > 0
                    IMG.SP(old_label).neighbors(neighbor_label) = num_neighbors-1;
                else
                    disp('Trying to remove a neighbor that was never added');
                end
            end
        end

        % update the neighbors' neighbors list. (not a typo!)
        if neighbor_labels(1)>0 && neighbor_labels(1)~=old_label
            IMG.SP(neighbor_labels(1)) = SP_update_neighbors_label_rem_check(IMG.SP(neighbor_labels(1)), IMG.label, x-1, y, old_label);
        end
        if neighbor_labels(2)>0 && neighbor_labels(2)~=old_label
            IMG.SP(neighbor_labels(2)) = SP_update_neighbors_label_rem_check(IMG.SP(neighbor_labels(2)), IMG.label, x+1, y, old_label);
        end
        if neighbor_labels(3)>0 && neighbor_labels(3)~=old_label
            IMG.SP(neighbor_labels(3)) = SP_update_neighbors_label_rem_check(IMG.SP(neighbor_labels(3)), IMG.label, x, y-1, old_label);
        end
        if neighbor_labels(4)>0 && neighbor_labels(4)~=old_label
            IMG.SP(neighbor_labels(4)) = SP_update_neighbors_label_rem_check(IMG.SP(neighbor_labels(4)), IMG.label, x, y+1, old_label);
        end
    end
end