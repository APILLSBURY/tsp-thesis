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
        if x>1
            llabel = IMG.label(x-1, y);
        else
            llabel = old_label;
        end
        if x<IMG.xdim
            rlabel = IMG.label(x+1, y);
        else
            rlabel = old_label;
        end
        if y>1
            ulabel = IMG.label(x, y-1);
        else
            ulabel = old_label;
        end
        if y<IMG.ydim
            dlabel = IMG.label(x, y+1);
        else
            dlabel = old_label;
        end

        if (llabel > 0 && llabel~=old_label)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), llabel);
        end
        if (rlabel > 0 && rlabel~=old_label && rlabel ~=llabel)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), rlabel);
        end
        if (ulabel > 0 && ulabel~=old_label && ulabel~=llabel && ulabel~=rlabel)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), ulabel);
        end
        if (dlabel > 0 && dlabel~=old_label && dlabel~=llabel && dlabel~=rlabel && dlabel~=ulabel)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), dlabel);
        end

        % update the neighbors' neighbors list. (not a typo!)
        if (llabel > 0 && llabel~=old_label)
            IMG.SP(llabel) = SP_update_neighbors_label_rem_check(IMG.SP(llabel), IMG.label, x-1, y, old_label);
        end
        if (rlabel > 0 && rlabel~=old_label)
            IMG.SP(rlabel) = SP_update_neighbors_label_rem_check(IMG.SP(rlabel), IMG.label, x+1, y, old_label);
        end
        if (ulabel > 0 && ulabel~=old_label)
            IMG.SP(ulabel) = SP_update_neighbors_label_rem_check(IMG.SP(ulabel), IMG.label, x, y-1, old_label);
        end
        if (dlabel > 0 && dlabel~=old_label)
            IMG.SP(dlabel) = SP_update_neighbors_label_rem_check(IMG.SP(dlabel), IMG.label, x, y+1, old_label);
        end
    end
end