% --------------------------------------------------------------------------
% -- update_neighbors_label_add
% --   Increments the corresponding neighbor count for neighbor_label
% --
% --   parameters:
% --     - neighbor_label : the label of the neighbor to increment
% --------------------------------------------------------------------------
function SP = SP_update_neighbors_label_add2(SP, neighbor_label)
    SP.neighbors(neighbor_label) = SP.neighbors(neighbor_label) + 1;
end