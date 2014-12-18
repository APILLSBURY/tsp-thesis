% --------------------------------------------------------------------------
% -- U_update_neighbor_list
% --   Updates super pixel neighbors used in merge
% --
% --   parameters:
% --     - neighbors : a boolean array indicating which labels are added
% --     - neighborsLL : a linked list of the neighbor indices
% --     - index : the index of the point to consider adding
% --------------------------------------------------------------------------
function neighbors = U_update_neighbor_list(IMG, neighbors, index)
    [x, y] = get_x_and_y_from_index(index, IMG.xdim);
    k = IMG.label(x, y);
    if k>0
        neighbors(k) = true;
    end
end