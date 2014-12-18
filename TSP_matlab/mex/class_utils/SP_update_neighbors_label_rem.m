% --------------------------------------------------------------------------
% -- update_neighbors_label_rem
% --   Decrements the corresponding neighbor count for neighbor_label
% --
% --   parameters:
% --     - neighbor_label : the label of the neighbor to decrement
% --------------------------------------------------------------------------
function SP = SP_update_neighbors_label_rem(SP, neighbor_label)
    if SP.neighbors(neighbor_label) > 0
        SP.neighbors(neighbor_label) = SP.neighbors(neighbor_label) - 1;
    else
        disp('Trying to remove a neighbor that was never added');
    end
end