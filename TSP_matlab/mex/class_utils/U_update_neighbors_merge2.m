% --------------------------------------------------------------------------
% -- U_update_neighbors_merge
% --   Updates the neighbor lists for merging two super pixels.
% -- All labels and SP_arr things should be updated *before* calling this.
% --
% --   parameters:
% --     - index : the index bordering the removed pixel
% --------------------------------------------------------------------------
function IMG = U_update_neighbors_merge2(IMG, new_label, old_label)
    IMG.SP(new_label).neighbors = IMG.SP(new_label).neighbors + IMG.SP(old_label).neighbors;
    
    %The new label can't have itself as its neighbor, and the old label is empty
    IMG.SP(new_label).neighbors(new_label) = 0;
    IMG.SP(new_label).neighbors(old_label) = 0;

    % now update all neighboring neighbor lists
    for neighbor_k=find(IMG.SP(new_label).neighbors)';
        if neighbor_k > 0
            IMG.SP(neighbor_k) = SP_update_neighbors_self(IMG.SP(neighbor_k), IMG.label);
        end
    end
end