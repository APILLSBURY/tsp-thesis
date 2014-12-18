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
        if x>1
            llabel = IMG.label(x-1, y);
        else
            llabel = cur_label;
        end
        if x<IMG.xdim
            rlabel = IMG.label(x+1, y);
        else
            rlabel = cur_label;
        end
        if y>1
            ulabel = IMG.label(x, y-1);
        else
            ulabel = cur_label;
        end
        if y<IMG.ydim
            dlabel = IMG.label(x, y+1);
        else
            dlabel = cur_label;
        end

        if (llabel > 0 && llabel~=cur_label)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), llabel);
        end
        if (rlabel > 0 && rlabel~=cur_label && rlabel ~=llabel)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), rlabel);
        end
        if (ulabel > 0 && ulabel~=cur_label && ulabel~=llabel && ulabel~=rlabel)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), ulabel);
        end
        if (dlabel > 0 && dlabel~=cur_label && dlabel~=llabel && dlabel~=rlabel && dlabel~=ulabel)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), dlabel);
        end

        % update the neighbors' neighbors list. (not a typo!)
        if (llabel > 0 && llabel~=cur_label)
            IMG.SP(llabel) = SP_update_neighbors_label_add_check(IMG.SP(llabel), IMG.label, x-1, y, cur_label);
        end
        if (rlabel > 0 && rlabel~=cur_label)
            IMG.SP(rlabel) = SP_update_neighbors_label_add_check(IMG.SP(rlabel), IMG.label, x+1, y, cur_label);
        end
        if (ulabel > 0 && ulabel~=cur_label)
            IMG.SP(ulabel) = SP_update_neighbors_label_add_check(IMG.SP(ulabel), IMG.label, x, y-1, cur_label);
        end
        if (dlabel > 0 && dlabel~=cur_label)
            IMG.SP(dlabel) = SP_update_neighbors_label_add_check(IMG.SP(dlabel), IMG.label, x, y+1, cur_label);
        end
    end
end