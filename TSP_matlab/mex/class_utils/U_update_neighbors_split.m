% --------------------------------------------------------------------------
% -- U_update_neighbors_split
% --   Updates the neighbor lists for merging two super pixels.
% -- All labels and SP_arr things should be updated *before* calling this.
% --
% --   parameters:
% --     - index : the index bordering the removed pixel
% --------------------------------------------------------------------------
function IMG = U_update_neighbors_split(IMG, label1, label2)
    IMG.SP(label1) = SP_update_neighbors_self(IMG.SP(label1), IMG.label);
    IMG.SP(label2) = SP_update_neighbors_self(IMG.SP(label2), IMG.label);

    for neighbor_k=1:length(IMG.SP(label1).neighbors)
        if (IMG.SP(label1).neighbors(neighbor_k) > 0 && neighbor_k ~= label2) ...
                || (IMG.SP(label2).neighbors(neighbor_k) > 0 && neighbor_k ~= label1)
            IMG.SP(neighbor_k) = SP_update_neighbors_self(IMG.SP(neighbor_k), IMG.label);
        end
    end
end
