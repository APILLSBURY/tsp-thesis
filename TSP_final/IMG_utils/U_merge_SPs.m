% --------------------------------------------------------------------------
% -- merge_with
% --   Merges this SP with another one and empties the other one. Assumes
% -- the labels have already been updated correctly. Fixes borders.
% --
% --   parameters:
% --     - other : A pointer to the other SP to merge with
% --     - label : A pointer to the label image
% --     - border_ptr : the linked list borders to fix
% --     - (xdim, ydim) : image dimensions
% --------------------------------------------------------------------------
function IMG_SP = U_merge_SPs(IMG_SP, IMG_label, index, index_other)
    % update the MIW stuff
    IMG_SP(index).N = IMG_SP(index).N + IMG_SP(index_other).N;
    IMG_SP(index).pixels = IMG_SP(index).pixels + IMG_SP(index_other).pixels;
    IMG_SP(index).app = NormalD_merge(IMG_SP(index).app, IMG_SP(index_other).app);
    IMG_SP(index).pos = NormalD_merge(IMG_SP(index).pos, IMG_SP(index_other).pos);
    
    % update the border pixels
    IMG_SP(index) = SP_fix_borders(IMG_SP(index), IMG_label);
    
    %update neighbors' neighbors lists
    found_neighbors = find(IMG_SP(index_other).neighbors);
    for nindex=1:length(found_neighbors)
        n = found_neighbors(nindex);
        if n>numel(IMG_SP)
            disp('wat is happening');
            disp(n);
            disp(numel(IMG_SP));
        elseif IMG_SP(n).neighbors(index)>0
            IMG_SP(n).neighbors = zeros(size(IMG_SP(n).neighbors));
            found_borders = find(IMG_SP(n).borders);
            for i=1:length(found_borders)
                IMG_SP(n) = SP_update_neighbors_add_self(IMG_SP(n), IMG_label, found_borders(i));
            end
        else
            IMG_SP(n).neighbors(index) = IMG_SP(n).neighbors(index_other);
            IMG_SP(n).neighbors(index_other) = 0;
        end
    end
    
    %update index and index_other's neighbor lists
    IMG_SP(index).neighbors = IMG_SP(index).neighbors + IMG_SP(index_other).neighbors;
    IMG_SP(index).neighbors(index) = 0;
    IMG_SP(index).neighbors(index_other) = 0;
    
    IMG_SP(index) = SP_calculate_log_probs(IMG_SP(index));
    IMG_SP(index_other) = SP_empty(IMG_SP(index_other));
end