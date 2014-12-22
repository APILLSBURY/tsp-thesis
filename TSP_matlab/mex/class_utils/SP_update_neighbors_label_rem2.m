% --------------------------------------------------------------------------
% -- update_neighbors_label_rem
% --   Decrements the corresponding neighbor count for neighbor_label
% --
% --   parameters:
% --     - neighbor_label : the label of the neighbor to decrement
% --------------------------------------------------------------------------
function SP = SP_update_neighbors_label_rem2(SP, neighbor_label)
    neighbor_val = SP.neighbors(neighbor_label);
    if neighbor_val > 0
        SP.neighbors(neighbor_label) = neighbor_val - 1;
    else
        disp('Trying to remove a neighbor that was never added');
    end
end