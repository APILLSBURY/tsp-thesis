% --------------------------------------------------------------------------
% -- move_local
% --   finds the optimal local joint move of labels and parameters. Chooses
% -- a super pixel at random and loops through and updates its borders.
% --------------------------------------------------------------------------
function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old, changed] = local_move_internal(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs)
    [xdim, ydim] = size(IMG_label);
    if (IMG_K>length(IMG_SP_changed))
        disp('Ran out of space!');
    end

    changed = false;

    % choose a random order of super pixels
    perm = zeros(1, IMG_max_SPs);
    perm(1, 1:IMG_K) = randperm(IMG_K);
    Nsp = IMG_K;
    % MM STEP
    for for_kindex=1:length(perm)
        for_k = perm(for_kindex);
        if for_k ~= 0
            % find a nonempty super pixel
            if ~(for_k > numel(IMG_SP) || isempty(IMG_SP(for_k).N) || IMG_SP(for_k).N == 0) && IMG_SP_changed(for_k)
                IMG_SP_changed(for_k) = false;

                % loop through borders
                found_borders = find(IMG_SP(for_k).borders)';
                border_length = length(found_borders);
                bfs_length = 4 * border_length;
                big_bfs_queue = zeros(bfs_length,1);
                big_bfs_index = 1;

                i = 1;
                while i<=border_length+bfs_length
                    if for_k > 0 && (for_k > numel(IMG_SP) || isempty(IMG_SP(for_k).N) || IMG_SP(for_k).N == 0)
                        break;
                    end
                    if i<=border_length
                        k = for_k;
                        index = found_borders(i);
                        
                        %make sure that it's still a border
                        if ~IMG_SP(k).borders(index)
                            break;
                        end
                        [x, y] = get_x_and_y_from_index(index, xdim);
                    else
                        if big_bfs_queue(i - border_length)>0
                            index = big_bfs_queue(i - border_length);
                        else
                            break;
                        end
                        [x, y] = get_x_and_y_from_index(index, xdim);
                        k = IMG_label(x, y);
                        if k>0
                            i = i+1;
                            continue;
                        end
                    end
                    i = i+1;

                    % temporarily remove this data point from the SP
                    max_prob = -inf;
                    max_k = -2;

                    % the current k has to be a possible choice
                    [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, index, IMG_label(x, y), false, max_prob, max_k);

                    if (~IMG_boundary_mask(x, y))
                        [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, index, 0, true, max_prob, max_k);
                    end

                    % find which k's we can move to
                    if (x>1 && IMG_label(x-1, y)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, index, IMG_label(x-1, y), true, max_prob, max_k);
                    end
                    if (y>1 && IMG_label(x, y-1)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, index, IMG_label(x, y-1), true, max_prob, max_k);
                    end
                    if (x<xdim && IMG_label(x+1, y)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, index, IMG_label(x+1, y), true, max_prob, max_k);
                    end
                    if (y<ydim && IMG_label(x, y+1)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, index, IMG_label(x, y+1), true, max_prob, max_k);
                    end


                    if (max_k==k && k>0)
                        IMG_SP(k).borders(index) = true;
                    elseif max_k~=k
                        changed = true;
                        top = y>1 && IMG_label(x,y-1)==k;
                        bottom = y<ydim && IMG_label(x,y+1)==k;
                        left = x>1 && IMG_label(x-1,y)==k;
                        right = x<xdim && IMG_label(x+1,y)==k;
                        if for_kindex <= Nsp && ((top && bottom && ~left && ~right) || (left && right && ~top && ~bottom)) && ...
                                k>0 && k<=numel(IMG_SP) && any(IMG_SP(k).neighbors)
                            %fprintf('doin a splitty thing with %d and %d\n', k, IMG_K+1);
                            IMG_K=IMG_K+1;
                            IMG_SP(IMG_K) = new_SP(IMG_new_pos, IMG_new_app, 0, [0, 0], IMG_N, IMG_max_SPs);
                            
                            cont1 = true;
                            cont2 = true;
                            fat1 = false;
                            fat2 = false;
                            
                            prev_x1 = x;
                            prev_x2 = x;
                            prev_y1 = y;
                            prev_y2 = y;
                            if left
                                x1 = x+1;
                                y1 = y;
                                x2 = x-1;
                                y2 = y;
                            else
                                x1 = x;
                                y1 = y+1;
                                x2 = x;
                                y2 = y-1;
                            end
                            sides = false(4, 2);
                            count = 0;
                            while cont1 || cont2
                                sides(1,1) = y1>1 && IMG_label(x1,y1-1)==k;
                                sides(2,1) = y1<ydim && IMG_label(x1,y1+1)==k;
                                sides(3,1) = x2>1 && IMG_label(x1-1,y1)==k;
                                sides(4,1) = x2<xdim && IMG_label(x1+1,y1)==k;
                                
                                sides(1,2) = y1>1 && IMG_label(x2,y2-1)==k;
                                sides(2,2) = y1<ydim && IMG_label(x2,y2+1)==k;
                                sides(3,2) = x2>1 && IMG_label(x2-1,y2)==k;
                                sides(4,2) = x2<xdim && IMG_label(x2+1,y2)==k;
                                if cont1
                                    count = count + 1;
                                    side_sum = sum(sides(:,1));
                                    if side_sum == 1
                                        cont1 = false;
                                        %disp('end1 is skinny');
                                    elseif side_sum == 2
                                        if sides(1,1) && prev_y1~=y1-1
                                            prev_y1 = y1;
                                            prev_x1 = x1;
                                            y1 = y1-1;
                                        elseif sides(2,1) && prev_y1~=y1+1
                                            prev_y1 = y1;
                                            prev_x1 = x1;
                                            y1 = y1+1;
                                        elseif sides(3,1) && prev_x1~=x1-1
                                            prev_y1 = y1;
                                            prev_x1 = x1;
                                            x1 = x1-1;
                                        else
                                            prev_y1 = y1;
                                            prev_x1 = x1;
                                            x1 = x1+1;
                                        end
                                    else
                                        cont1 = false;
                                        fat1 = true;
                                        %disp('end1 is fat');
                                    end
                                end
                                if cont2
                                    count = count + 1;
                                    side_sum = sum(sides(:,2));
                                    if side_sum == 1
                                        cont2 = false;
                                        %disp('end2 is skinny');
                                    elseif side_sum == 2
                                        if sides(1,2) && prev_y2~=y2-1
                                            prev_y2 = y2;
                                            prev_x2 = x2;
                                            y2 = y2-1;
                                        elseif sides(2,2) && prev_y2~=y2+1
                                            prev_y2 = y2;
                                            prev_x2 = x2;
                                            y2 = y2+1;
                                        elseif sides(3,2) && prev_x2~=x2-1
                                            prev_y2 = y2;
                                            prev_x2 = x2;
                                            x2 = x2-1;
                                        else
                                            prev_y2 = y2;
                                            prev_x2 = x2;
                                            x2 = x2+1;
                                        end
                                    else
                                        fat2 = true;
                                        cont2 = false;
                                        %disp('end2 is fat');
                                    end
                                end
                            end
                            
                            if count > 10
                                pos_mean = NormalD_calc_mean(IMG_SP(k).pos);
                                dist1 = sqrt((pos_mean(1)-x1)^2 + (pos_mean(2)-y1)^2);
                                dist2 = sqrt((pos_mean(1)-x2)^2 + (pos_mean(2)-y2)^2);
                                starting = [prev_x1 prev_x2; prev_y1 prev_y2; x1 x2; y1 y2];
                                if fat1 && (~fat2 || dist1<dist2)
                                    starting_index = 1;
                                else %either fat2 or both are skinny (which probably won't happen much)
                                    starting_index = 2;
                                end

                                fprintf('bfs_splitting on %d and %d\n', k, IMG_K);
                                [IMG_SP, IMG_label] = bfs_split(IMG_SP, IMG_label, IMG_N, IMG_data, IMG_boundary_mask, k, IMG_K, starting, starting_index, ~(fat1 && fat2));

                                %if the new TSP is larger, switch the new and old
                                if IMG_SP(IMG_K).N > IMG_SP(k).N
                                    fprintf('doin a switchy thing now! IMG_SP(%d).N is %d, IMG_SP(%d).N is %d\n', IMG_K, IMG_SP(IMG_K).N, k, IMG_SP(k).N);
                                    if IMG_SP(k).N == 0
                                        disp('um excuse me why is it 0');
                                    end
                                    
                                    %switch the pixel labels
                                    IMG_label(IMG_label==k) = IMG_K+1;
                                    IMG_label(IMG_label==IMG_K) = k;
                                    IMG_label(IMG_label==IMG_K+1) = IMG_K;
                                    
                                    %switch the other info by merging to a
                                    %temporary SP
                                    if IMG_K+1>IMG_max_SPs
                                        disp('Ran out of space!');
                                    end
                                    IMG_SP(IMG_K+1) = new_SP(IMG_new_pos, IMG_new_app, IMG_max_UID, [0, 0], IMG_N, IMG_max_SPs);
                                    IMG_SP = U_merge_SPs(IMG_SP, IMG_label, IMG_K+1, k);
                                    IMG_SP = U_merge_SPs(IMG_SP, IMG_label, k, IMG_K);
                                    IMG_SP = U_merge_SPs(IMG_SP, IMG_label, IMG_K, IMG_K+1);

                                    %delete the temp SP
                                    IMG_SP(IMG_K+1) = [];
                                end

                                %move every pixel in the new SP to a different SP
                                while IMG_SP(IMG_K).N > 0
                                    temp_found_borders = find(IMG_SP(IMG_K).borders);
                                    for temp_border_index=1:length(temp_found_borders)
                                        curr_index = temp_found_borders(temp_border_index);
                                        [temp_x, temp_y] = get_x_and_y_from_index(curr_index, xdim);
                                        max_k = -2;
                                        max_prob = -inf;
                                        % find which k's we can move to
                                        if (x>1 && IMG_label(temp_x-1, temp_y)~=IMG_K && IMG_label(temp_x-1, temp_y)>0)
                                            [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, curr_index, IMG_label(temp_x-1, temp_y), true, max_prob, max_k);
                                        end
                                        if (y>1 && IMG_label(temp_x, temp_y-1)~=IMG_K && IMG_label(temp_x, temp_y-1)>0)
                                            [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, curr_index, IMG_label(temp_x, temp_y-1), true, max_prob, max_k);
                                        end
                                        if (x<xdim && IMG_label(temp_x+1, temp_y)~=IMG_K && IMG_label(temp_x+1, temp_y)>0)
                                            [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, curr_index, IMG_label(temp_x+1, temp_y), true, max_prob, max_k);
                                        end
                                        if (y<ydim && IMG_label(temp_x, temp_y+1)~=IMG_K && IMG_label(temp_x, temp_y+1)>0)
                                            [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, curr_index, IMG_label(temp_x, temp_y+1), true, max_prob, max_k);
                                        end

                                        %check to make sure we're bordering another SP, not just bordering the edge of the image
                                        if max_k~=-2 && max_prob ~= -inf
                                            IMG_label(temp_x, temp_y) = max_k;

                                            % update the neighbors lists
                                            IMG_SP = U_update_neighbors_rem(IMG_label, IMG_SP, IMG_K, curr_index);
                                            IMG_SP = U_update_neighbors_add(IMG_label, IMG_SP, curr_index);

                                            % set the correct SP_changed variables of all neighbors
                                            IMG_SP_changed(IMG_SP(IMG_K).neighbors > 0) = true;
                                            IMG_SP_changed(IMG_SP(max_k).neighbors > 0) = true;

                                            % update all border lists for neighbors
                                            IMG_SP = U_update_border_changed(IMG_label, IMG_SP, curr_index);

                                            % add this point to the maximum SP
                                            IMG_SP(IMG_K) = SP_rem_pixel(IMG_SP(IMG_K), IMG_data, curr_index, IMG_boundary_mask(temp_x, temp_y));
                                            IMG_SP(max_k) = SP_add_pixel(IMG_SP(max_k), IMG_data, curr_index, U_check_border_pix(IMG_label, curr_index), IMG_boundary_mask(temp_x, temp_y));
                                        end
                                    end
                                end
                            end
                            
                            % remove the new SP now that it's empty
                            IMG_SP(IMG_K) = [];
                            perm(IMG_K) = k;
                            IMG_K = IMG_K-1;
                        
                        % check the topology for this pixel
                        elseif check_topology(IMG_label, index, IMG_T4Table)
                            %fprintf('making a local move from %d to %d\n', k, max_k);
                            % update the labels... it moves from k->max_k
                            IMG_label(x, y) = max_k;
                            if (max_k>IMG_K)
                                % creating a new one
                                disp('Creating a new TSP in local_move_internal');
                                if (IMG_K>=IMG_N)
                                    disp('Ran out of space!');
                                end
                                IMG_SP(IMG_K+1) = new_SP(IMG_new_pos, IMG_new_app, IMG_max_UID, [0, 0], IMG_N, IMG_max_SPs);
                                IMG_max_UID = IMG_max_UID + 1;
                                IMG_K = IMG_K + 1;
                            end

                            % update the neighbors lists
                            IMG_SP = U_update_neighbors_rem(IMG_label, IMG_SP, k, index);
                            IMG_SP = U_update_neighbors_add(IMG_label, IMG_SP, index);

                            % set the correct SP_changed variables of all neighbors
                            if (k<1)
                                IMG_SP_changed(max_k) = true;
                            else
                                IMG_SP_changed(IMG_SP(k).neighbors > 0) = true;
                                if (max_k>0)
                                    IMG_SP_changed(IMG_SP(max_k).neighbors > 0) = true;
                                end
                            end


                            % update all border lists for neighbors
                            IMG_SP = U_update_border_changed(IMG_label, IMG_SP, index);
                            if (k>0)
                                IMG_SP(k) = SP_rem_pixel(IMG_SP(k), IMG_data, index, IMG_boundary_mask(x, y));
                            end

                            if k>0 && k <= numel(IMG_SP) && (isempty(IMG_SP(k).N) || IMG_SP(k).N == 0)
                                if (IMG_SP_old(k))
                                    IMG_alive_dead_changed = true;
                                else
                                    fprintf('killing TSP %d in local_move_internal\n', k);
                                    IMG_SP(k) = SP_empty(IMG_SP(k));
                                    IMG_SP_changed(k) = false;
                                end
                            end

                            % add this point to the maximum SP
                            if max_k>0
                                IMG_SP(max_k) = SP_add_pixel(IMG_SP(max_k), IMG_data, index, U_check_border_pix(IMG_label, index), IMG_boundary_mask(x, y));
                            end
                        end
                    end
                    if (x>1 && IMG_label(x-1, y)<1 && big_bfs_index <= bfs_length)
                        big_bfs_queue(big_bfs_index) = index-1;
                        big_bfs_index = big_bfs_index + 1;
                    end
                    if (y>1 && IMG_label(x, y-1)<1 && big_bfs_index <= bfs_length)
                        big_bfs_queue(big_bfs_index) = index-xdim;
                        big_bfs_index = big_bfs_index + 1;
                    end
                    if (x<xdim && IMG_label(x+1, y)<1 && big_bfs_index <= bfs_length)
                        big_bfs_queue(big_bfs_index) = index+1;
                        big_bfs_index = big_bfs_index + 1;
                    end
                    if (y<ydim && IMG_label(x, y+1)<1 && big_bfs_index <= bfs_length)
                        big_bfs_queue(big_bfs_index) = index+xdim;
                        big_bfs_index = big_bfs_index + 1;
                    end
                end
            end
        end
    end
    if ~changed
        IMG_SP_changed(1:IMG_K) = false;
    end
end


% --------------------------------------------------------------------------
% -- move_local_calc_neigbor
% --   calculates the probability of assigning the pixel at index to the
% -- cluster of nindex (neighbor index). Updates the max_prob and max_k if
% -- the new_app value is greater than the old one
% --
% --   parameters:
% --     - index : the new point to add
% --     - new_k : the neighbor to add this point to
% --     - max_prob : the maximum probability of all neighbors
% --     - max_k : the super pixel index of the maximum probability
% --------------------------------------------------------------------------
function [max_prob, max_k] = move_local_calc_delta(IMG_label, IMG_boundary_mask, IMG_dummy_log_prob, IMG_SP, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, index, new_k, add, max_prob, max_k)
    [xdim, ~] = size(IMG_label);
    prob = 0;
    [x, y] = get_x_and_y_from_index(index, xdim);

    if (new_k<1)
        if IMG_boundary_mask(x,y)
            prob = -inf;
        else
            prob = IMG_dummy_log_prob;
        end
    else
        if (new_k > numel(IMG_SP) || isempty(IMG_SP(new_k).N) || IMG_SP(new_k).N == 0)
            temp_SP = IMG_new_SP;
            is_old = false;
        else
            temp_SP = IMG_SP(new_k);
            is_old = IMG_SP_old(new_k);
        end

        SP_size = temp_SP.N;
        if (add)
            prob = prob + SP_log_likelihood_test_point(temp_SP, IMG_data(index,:), IMG_boundary_mask(x, y));
            prob = prob - temp_SP.log_likelihood;
            if SP_size>0
                prob = prob + (model_order_params.area - (SP_size+1)/2)*(SP_size+1)/model_order_params.area_var ...
                    - (model_order_params.area - SP_size/2)*SP_size/model_order_params.area_var;
            else % size == 0
                prob = prob + (model_order_params.area - 1/2)/model_order_params.area_var;
                if is_old
                    prob = prob + model_order_params.is_old_const;
                else
                    prob = prob + model_order_params.is_new_const;
                end
            end
        else
            prob = prob + temp_SP.log_likelihood;
            prob = prob - SP_log_likelihood_test_point_rem(temp_SP, IMG_data(index,:), IMG_boundary_mask(x, y));

            if SP_size>1
                prob = prob + (model_order_params.area - SP_size/2)*SP_size/model_order_params.area_var ...
                    - (model_order_params.area - (SP_size-1)/2)*(SP_size-1)/model_order_params.area_var;
            elseif SP_size==1
                prob = prob + (model_order_params.area - 1/2)/model_order_params.area_var;
                if is_old
                    prob = prob + model_order_params.is_old_const;
                else
                    prob = prob + model_order_params.is_new_const;
                end
            end % else if size==0 then no change
        end
    end

    if (prob>max_prob+1e-10 || max_k==-2)
        max_prob = prob;
        max_k = new_k;
    end
end

function top_OK = check_topology(IMG_label, index, IMG_T4Table)
    [xdim, ydim] = size(IMG_label);
    [x, y] = get_x_and_y_from_index(index, xdim);
    values = ones(4,1);
    if x>1
        values(1) = topology_number_label(IMG_label, x, y, IMG_label(x-1, y), IMG_T4Table);
    end
    if y>1
        values(2) = topology_number_label(IMG_label, x, y, IMG_label(x, y-1), IMG_T4Table);
    end
    if x<xdim
        values(3) = topology_number_label(IMG_label, x, y, IMG_label(x+1, y), IMG_T4Table);
    end
    if y<ydim
        values(4) = topology_number_label(IMG_label, x, y, IMG_label(x, y+1), IMG_T4Table);
    end
    top_OK = all(values==1);
end

function value = topology_number_label(IMG_label, x, y, check_label, IMG_T4Table)
    [xdim, ydim] = size(IMG_label);
    neighborhood = false(8,1);
    neighborhood(2) = y>1 && IMG_label(x, y-1)==check_label;
    neighborhood(4) = x<xdim && IMG_label(x+1, y)==check_label;
    neighborhood(6) = y<ydim && IMG_label(x, y+1)==check_label;
    neighborhood(8) = x>1 && IMG_label(x-1, y)==check_label;
    neighborhood(1) = (neighborhood(2) || neighborhood(8)) && (x>1 && y>1 && IMG_label(x-1, y-1)==check_label);
    neighborhood(3) = (neighborhood(2) || neighborhood(4)) && (x<xdim && y>1 && IMG_label(x+1, y-1)==check_label);
    neighborhood(5) = (neighborhood(4) || neighborhood(6)) && (x<xdim && y<ydim && IMG_label(x+1, y+1)==check_label);
    neighborhood(7) = (neighborhood(6) || neighborhood(8)) && (x>1 && y<ydim && IMG_label(x-1, y+1)==check_label);

    index = 0;
    for i=8:-1:1
        index = index*2;
        if neighborhood(i)
            index = index+1;
        end
    end
    value = IMG_T4Table(index+1);
end

function [IMG_SP, IMG_label] = bfs_split(IMG_SP, IMG_label, IMG_N, IMG_data, IMG_boundary_mask, k, IMG_K, starting, starting_index, switched);
    %set the prev values to temp and the x/y values to prev
    temp_x = starting(1, starting_index);
    temp_y = starting(2, starting_index);
    prev_x = starting(3, starting_index);
    prev_y = starting(4, starting_index);
    
    bfs_write = 1;
    bfs_read = 1;
    bfs = zeros(IMG_N, 2);
    [xdim, ydim] = size(IMG_label);
    visited = false(xdim, ydim);
    max_pixels = min(max(100, IMG_SP(k).N*0.02), IMG_SP(k).N/2);
    
    if starting_index == 1
        other_side = [starting(1, 2), starting(2, 2)];
    else
        other_side = [starting(1, 1), starting(2, 1)]; 
    end
    
    %change prev_x, prev_y temporarily and then change it back so we can
    %check for loops
    changed_back = false;
    IMG_label(prev_x, prev_y) = IMG_K;
    prev_index = get_index_from_x_and_y(prev_x, prev_y, xdim);
    IMG_SP = U_update_neighbors_rem(IMG_label, IMG_SP, k, prev_index);
    IMG_SP = U_update_neighbors_add(IMG_label, IMG_SP, prev_index);
    IMG_SP(k) = SP_rem_pixel(IMG_SP(k), IMG_data, prev_index, IMG_boundary_mask(prev_x, prev_y));
    IMG_SP(IMG_K) = SP_add_pixel(IMG_SP(IMG_K), IMG_data, prev_index, U_check_border_pix(IMG_label, prev_index), IMG_boundary_mask(prev_x, prev_y));
    IMG_SP = U_update_border_changed(IMG_label, IMG_SP, prev_index);
%     visited(prev_x, prev_y) = true;

    bfs(bfs_write, :) = [temp_x, temp_y];
    visited(temp_x, temp_y) = true;
    bfs_write = bfs_write+1;

    % switch all connected pixels from k to IMG_K
    while bfs_write > bfs_read
        temp_coor = bfs(bfs_read, :);
        temp_x = temp_coor(1);
        temp_y = temp_coor(2);
        temp_index = get_index_from_x_and_y(temp_x, temp_y, xdim);

        IMG_label(temp_x, temp_y) = IMG_K;
        IMG_SP = U_update_neighbors_rem(IMG_label, IMG_SP, k, temp_index);
        IMG_SP = U_update_neighbors_add(IMG_label, IMG_SP, temp_index);
        IMG_SP(k) = SP_rem_pixel(IMG_SP(k), IMG_data, temp_index, IMG_boundary_mask(temp_x, temp_y));
        IMG_SP(IMG_K) = SP_add_pixel(IMG_SP(IMG_K), IMG_data, temp_index, U_check_border_pix(IMG_label, temp_index), IMG_boundary_mask(temp_x, temp_y));
        IMG_SP = U_update_border_changed(IMG_label, IMG_SP, temp_index);

        if (temp_x>1 && IMG_label(temp_x-1, temp_y)==k) && ~visited(temp_x-1, temp_y)
            visited(temp_x-1, temp_y) = true;
            bfs(bfs_write,:) = [temp_x-1, temp_y];
            bfs_write = bfs_write + 1;
        end
        if (temp_y>1 && IMG_label(temp_x, temp_y-1)==k) && ~visited(temp_x, temp_y-1)
            visited(temp_x, temp_y-1) = true;
            bfs(bfs_write,:) = [temp_x, temp_y-1];
            bfs_write = bfs_write + 1;
        end
        if (temp_x<xdim && IMG_label(temp_x+1, temp_y)==k) && ~visited(temp_x+1, temp_y)
            visited(temp_x+1, temp_y) = true;
            bfs(bfs_write,:) = [temp_x+1, temp_y];
            bfs_write = bfs_write + 1;
        end
        if (temp_y<ydim && IMG_label(temp_x, temp_y+1)==k) && ~visited(temp_x, temp_y+1)
            visited(temp_x, temp_y+1) = true;
            bfs(bfs_write,:) = [temp_x, temp_y+1];
            bfs_write = bfs_write + 1;
        end
        bfs_read = bfs_read+1;
        
        if ~changed_back
            changed_back = true;

            IMG_label(prev_x, prev_y) = k;
            IMG_SP = U_update_neighbors_rem(IMG_label, IMG_SP, IMG_K, prev_index);
            IMG_SP = U_update_neighbors_add(IMG_label, IMG_SP, prev_index);
            IMG_SP(IMG_K) = SP_rem_pixel(IMG_SP(IMG_K), IMG_data, prev_index, IMG_boundary_mask(prev_x, prev_y));
            IMG_SP(k) = SP_add_pixel(IMG_SP(k), IMG_data, prev_index, U_check_border_pix(IMG_label, prev_index), IMG_boundary_mask(prev_x, prev_y));
            IMG_SP = U_update_border_changed(IMG_label, IMG_SP, prev_index);
        end
        
        if visited(prev_x, prev_y)
            disp('we found a loop!');
            
            %reset all the pixels we changed
            IMG_label(visited) = k;
            IMG_SP = U_merge_SPs(IMG_SP, IMG_label, k, IMG_K);
            
            %block off the far side by setting it equal to IMG_K
            temp_index = get_index_from_x_and_y(other_side(1), other_side(2), xdim);
            IMG_label(other_side(1), other_side(2)) = IMG_K;
            IMG_SP = U_update_neighbors_rem(IMG_label, IMG_SP, k, temp_index);
            IMG_SP = U_update_neighbors_add(IMG_label, IMG_SP, temp_index);
            IMG_SP(k) = SP_rem_pixel(IMG_SP(k), IMG_data, temp_index, IMG_boundary_mask(other_side(1), other_side(2)));
            IMG_SP(IMG_K) = SP_add_pixel(IMG_SP(IMG_K), IMG_data, temp_index, U_check_border_pix(IMG_label, temp_index), IMG_boundary_mask(other_side(1), other_side(2)));
            IMG_SP = U_update_border_changed(IMG_label, IMG_SP, temp_index);
            
            %call bfs_split on just the skinny portion
            [IMG_SP, IMG_label] = bfs_split(IMG_SP, IMG_label, IMG_N, IMG_data, IMG_boundary_mask, k, IMG_K, starting, starting_index, true);
            break;
        end
        
        if ~switched && bfs_write > max_pixels
            fprintf('switching the switch on %d and %d\n', k, IMG_K);
            IMG_label(visited) = k;
            IMG_SP = U_merge_SPs(IMG_SP, IMG_label, k, IMG_K);
            if starting_index==1
                starting_index=2;
            else
                starting_index=1;
            end
            [IMG_SP, IMG_label] = bfs_split(IMG_SP, IMG_label, IMG_N, IMG_data, IMG_boundary_mask, k, IMG_K, starting, starting_index, true);
            break;
        end
    end
end