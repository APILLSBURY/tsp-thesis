% --------------------------------------------------------------------------
% -- U_update_neighbors_add
% --   Updates all neighbor lists when a pixel is "added". A previous
% -- call to update_neighbors_rem should be completed right before this one.
% -- The neighboring label should be changed before calling this function.
% --
% --   parameters:
% --     - index : the index bordering the removed pixel
% --------------------------------------------------------------------------
function IMG = U_update_neighbors_add(IMG, index)
    [x, y] = get_x_and_y_from_index(index, IMG.xdim);
    cur_label = IMG.label(x, y);
    if (cur_label>0)
        [x, y] = get_x_and_y_from_index(index, IMG.xdim);
        neighbor_labels = zeros(1, 4);
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
            if neighbor_label > 0 && neighbor_label ~= cur_label
                IMG.SP(cur_label).neighbors(neighbor_label) = IMG.SP(cur_label).neighbors(neighbor_label)+1;
            end
        end

        % update the neighbors' neighbors lists. (not a typo!)
        if neighbor_labels(1)>0 && neighbor_labels(1)~=cur_label
            IMG.SP(neighbor_labels(1)) = SP_update_neighbors_label_add_check(IMG.SP(neighbor_labels(1)), IMG.label, x-1, y, cur_label);
        end
        if neighbor_labels(2)>0 && neighbor_labels(2)~=cur_label
            IMG.SP(neighbor_labels(2)) = SP_update_neighbors_label_add_check(IMG.SP(neighbor_labels(2)), IMG.label, x+1, y, cur_label);
        end
        if neighbor_labels(3)>0 && neighbor_labels(3)~=cur_label
            IMG.SP(neighbor_labels(3)) = SP_update_neighbors_label_add_check(IMG.SP(neighbor_labels(3)), IMG.label, x, y-1, cur_label);
        end
        if neighbor_labels(4)>0 && neighbor_labels(4)~=cur_label
            IMG.SP(neighbor_labels(4)) = SP_update_neighbors_label_add_check(IMG.SP(neighbor_labels(4)), IMG.label, x, y+1, cur_label);
        end
    end
end