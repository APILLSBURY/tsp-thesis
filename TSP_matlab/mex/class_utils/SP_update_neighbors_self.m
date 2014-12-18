% --------------------------------------------------------------------------
% -- update_neighbors_self
% --   Updates the neighbor lists and counts by looking at all borders.
% -- Empties previous list. The borders list must be correct!
% --
% --   parameters:
% --     - label : the label image
% --     - (xdim,ydim) : dimensions of image
% --------------------------------------------------------------------------
function SP = SP_update_neighbors_self(SP, label)
    for index=1:length(SP.borders)
        if SP.borders(index)
            SP = SP_update_neighbors_add_self(SP, label, index);
        end
    end
end