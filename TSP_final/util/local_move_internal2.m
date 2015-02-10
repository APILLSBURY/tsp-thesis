% --------------------------------------------------------------------------
% -- move_local
% --   finds the optimal local joint move of labels and parameters. Chooses
% -- a super pixel at random and loops through and updates its borders.
% --------------------------------------------------------------------------
function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old, changed] = local_move_internal(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs)
    [xdim, ydim] = size(IMG_label);
    Nsp = numel(IMG_SP);
    if (IMG_K>length(IMG_SP_changed))
        disp('Ran out of space!');
    end

    changed = false;

    % temporary neighborhood
    neighborhood = false(9,1);

    % choose a random order of super pixels
    perm = randperm(IMG_K);
    
    % MM STEP
    for for_kindex=1:length(perm)
        for_k = perm(for_kindex);
        % find a nonempty super pixel
        if ~(for_k > numel(IMG_SP) || isempty(IMG_SP(for_k).N) || IMG_SP(for_k).N == 0) && IMG_SP_changed(for_k)
            IMG_SP_changed(for_k) = false;
            
            % loop through borders
            found_borders = find(IMG_SP(for_k).borders)';
            border_length = length(found_borders);
            bfs_length = 4 * border_length;
            bfs_queue = zeros(bfs_length,1);
            bfs_index = 1;
            
            i = 1;
            while i<=border_length+bfs_length
                if for_k > 0 && (for_k > numel(IMG_SP) || isempty(IMG_SP(for_k).N) || IMG_SP(for_k).N == 0)
                    break;
                end
                if i<=border_length
                    k = for_k;
                    index = found_borders(i);
                    [x, y] = get_x_and_y_from_index(index, xdim);
                else
                    if bfs_queue(i - border_length)>0
                        index = bfs_queue(i - border_length);
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
                    % check the topology for this pixel
                    if ~check_topology(IMG_label, index, neighborhood, IMG_T4Table)
                        if k>0 && k<=numel(IMG_SP) && any(IMG_SP(k).neighbors)
                            [IMG_K, IMG_label, IMG_SP, ~, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = move_split_SP(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_alive_dead_changed, IMG_N, IMG_new_SP, true(size(IMG_SP_changed)), Nsp, k);
                        end
                    else
                        fprintf('making a local move from %d to %d\n', k, max_k);
                        % update the labels... it moves from k->max_k
                        IMG_label(x, y) = max_k;
                        if (max_k>IMG_K)
                            % creating a new one
                            disp('Creating a new TSP in local_move_internal');
                            if (IMG_K>=IMG_N)
                                disp('Ran out of space!');
                            end
                            IMG_SP(IMG_K+1) = new_SP(IMG_new_pos, IMG_new_app, IMG_max_UID, [0, 0], IMG_N);
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
                                disp('killing a TSP in local_move_internal');
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
                if (x>1 && IMG_label(x-1, y)<1 && bfs_index <= bfs_length)
                    bfs_queue(bfs_index) = index-1;
                    bfs_index = bfs_index + 1;
                end
                if (y>1 && IMG_label(x, y-1)<1 && bfs_index <= bfs_length)
                    bfs_queue(bfs_index) = index-xdim;
                    bfs_index = bfs_index + 1;
                end
                if (x<xdim && IMG_label(x+1, y)<1 && bfs_index <= bfs_length)
                    bfs_queue(bfs_index) = index+1;
                    bfs_index = bfs_index + 1;
                end
                if (y<ydim && IMG_label(x, y+1)<1 && bfs_index <= bfs_length)
                    bfs_queue(bfs_index) = index+xdim;
                    bfs_index = bfs_index + 1;
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

function top_ok = check_topology(IMG_label, index, neighborhood, IMG_T4Table)
    [xdim, ydim] = size(IMG_label);
    [x, y] = get_x_and_y_from_index(index, xdim);
    top_ok = ( x<=1 || topology_number_label(IMG_label, x, y, IMG_label(x-1, y), neighborhood, IMG_T4Table)==1 ) && ...
            ( y<=1 || topology_number_label(IMG_label, x, y, IMG_label(x, y-1), neighborhood, IMG_T4Table)==1 ) && ...
            ( x>=xdim || topology_number_label(IMG_label, x, y, IMG_label(x+1, y), neighborhood, IMG_T4Table)==1 ) && ...
            ( y>=ydim || topology_number_label(IMG_label, x, y, IMG_label(x, y+1), neighborhood, IMG_T4Table)==1 );
end

function value = topology_number_label(IMG_label, x, y, check_label, neighborhood, IMG_T4Table)
    [xdim, ydim] = size(IMG_label);
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