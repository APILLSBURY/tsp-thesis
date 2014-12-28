% --------------------------------------------------------------------------
% -- move_local
% --   finds the optimal local joint move of labels and parameters. Chooses
% -- a super pixel at random and loops through and updates its borders.
% --------------------------------------------------------------------------
function IMG = local_move_internal(IMG)

    if (IMG.K>IMG.max_SPs)
        disp('Ran out of space!');
    end

    IMG.changed = false;

    % temporary neighborhood
    neighborhood = false(9,1);

    % choose a random order of super pixels
    permx = randperm(IMG.xdim-2*IMG.w) + IMG.w;
    permy = randperm(IMG.ydim-2*IMG.w) + IMG.w;
    
    % First check all the pixels on the actual image
    for x=permx
        for y=permy
            k = IMG.label(x, y);
            index = get_index_from_x_and_y(x, y, IMG.xdim);
            %%check to see if we're on a border
            if k>0 && IMG.SP(k).borders(index)
                % check the topology for this pixel
                if check_topology(IMG, index, neighborhood)
                    % temporarily remove this data point from the SP
                    max_prob = -inf;
                    max_k = -2;

                    % the current k has to be a possible choice
                    [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x, y), false, max_prob, max_k);

                    if (~IMG.boundary_mask(x, y))
                        [max_prob, max_k] = move_local_calc_delta(IMG, index, 0, true, max_prob, max_k);
                    end

                    % find which k's we can move to
                    if (x>1 && IMG.label(x-1, y)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x-1, y), true, max_prob, max_k);
                    end
                    if (y>1 && IMG.label(x, y-1)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x, y-1), true, max_prob, max_k);
                    end
                    if (x<IMG.xdim && IMG.label(x+1, y)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x+1, y), true, max_prob, max_k);
                    end
                    if (y<IMG.ydim && IMG.label(x, y+1)~=k)
                        [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x, y+1), true, max_prob, max_k);
                    end


                    if max_k~=k
                        IMG.changed = true;
                        % update the labels... it moves from k->max_k
                        IMG.label(x, y) = max_k;
                        if (max_k>IMG.K)
                            disp('Creating a new SP in local_move_internal');
                            % creating a new one
                            if (IMG.K>=IMG.max_SPs)
                                disp('Ran out of space!');
                            else
                                IMG.SP(IMG.K+1) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, [0, 0], IMG.N, IMG.max_SPs);
                                IMG.max_UID = IMG.max_UID + 1;
                                IMG.K = IMG.K + 1;
                                max_k = 0; % do this so the program breaks if this happens
                            end
                        end

                        
                        %remove the pixel from k and add it to max_k
                        if (k>0)
                            IMG.SP(k) = SP_rem_pixel(IMG.SP(k), IMG.data, index, IMG.boundary_mask(x, y));
                            IMG = U_update_neighbors_rem(IMG, k, index);
                        end
                        if max_k>0
                            IMG.SP(max_k) = SP_add_pixel(IMG.SP(max_k), IMG.data, index, U_check_border_pix(IMG, index), IMG.boundary_mask(x, y));
                            IMG = U_update_neighbors_add(IMG, index);
                        end
                        

                        % set the correct SP_changed variables of all neighbors
                        if (k<1)
                            IMG.SP_changed(max_k) = true;
                        else
                            IMG.SP_changed(IMG.SP(k).neighbors > 0) = true;
                            if (max_k>0)
                                IMG.SP_changed(IMG.SP(max_k).neighbors > 0) = true;
                            end
                        end

                        % update border lists for neighboring pixels
                        IMG = U_update_border_changed(IMG, index);

                        if k>0 && SP_is_empty(IMG, k)
                            if (IMG.SP_old(k))
                                IMG.alive_dead_changed = true;
                            else
                                IMG.SP_changed(k) = false;
                            end
                        end
                    end
                end
            end
        end
    end
        
%             % find a nonempty super pixel
%                 IMG.SP_changed(k) = false;
%                 % loop through borders
%                 border_length = length(found_borders);
%                 bfs_length = 5 * border_length;
%                 bfs_queue = zeros(bfs_length,1);
%                 bfs_index = 1;
% 
%                 i = 1;
%                 while (i<=border_length+bfs_length)
%                     if SP_is_empty(IMG, for_k)
%                         break;
%                     end
%                     if i<=border_length
%                         k = for_k;
%                         index = found_borders(i);
%                         [x, y] = get_x_and_y_from_index(index, IMG.xdim);
%                     else
%                         if bfs_queue(i - border_length) ~= 0
%                             index = bfs_queue(i - border_length);
%                         else
%                             break;
%                         end
%                         [x, y] = get_x_and_y_from_index(index, IMG.xdim);
%                         k = IMG.label(x, y);
%                         if k>0
%                             i = i+1;
%                             continue;
%                         end
%                     end
%                     i = i+1;
%                     
%                     do your thing
%                     
%                     
%                         if (x>1 && IMG.label(x-1, y)<1 && bfs_index <= bfs_length)
%                             bfs_queue(bfs_index) = index-1;
%                             bfs_index = bfs_index + 1;
%                         end
%                         if (y>1 && IMG.label(x, y-1)<1 && bfs_index <= bfs_length)
%                             bfs_queue(bfs_index) = index-IMG.xdim;
%                             bfs_index = bfs_index + 1;
%                         end
%                         if (x<IMG.xdim && IMG.label(x+1, y)<1 && bfs_index <= bfs_length)
%                             bfs_queue(bfs_index) = index+1;
%                             bfs_index = bfs_index + 1;
%                         end
%                         if (y<IMG.ydim && IMG.label(x, y+1)<1 && bfs_index <= bfs_length)
%                             bfs_queue(bfs_index) = index+IMG.xdim;
%                             bfs_index = bfs_index + 1;
%                         end
    if (~IMG.changed)
        IMG.SP_changed(1:IMG.K) = false;
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
function [max_prob, max_k] = move_local_calc_delta(IMG, index, new_k, add, max_prob, max_k)
    prob = 0;
    [x, y] = get_x_and_y_from_index(index, IMG.xdim);

    if (new_k<1)
        if IMG.boundary_mask(x,y)
            prob = -inf;
        else
            prob = IMG.dummy_log_prob;
        end
    else
        if SP_is_empty(IMG, new_k) % new super pixel
            temp_SP = create_SP_new(IMG);
            is_old = false;
        else
            temp_SP = IMG.SP(new_k);
            is_old = IMG.SP_old(new_k);
        end

        if (add)
            prob = prob + SP_log_likelihood_test_point(temp_SP, IMG.data(index,:), IMG.boundary_mask(x, y));
            prob = prob - temp_SP.log_likelihood;

            prob = prob + U_calc_model_order(IMG, temp_SP.N+1, is_old);
            prob = prob - U_calc_model_order(IMG, temp_SP.N, is_old);
        else
            prob = prob + temp_SP.log_likelihood;
            prob = prob - SP_log_likelihood_test_point_rem(temp_SP, IMG.data(index,:), IMG.boundary_mask(x, y));

            prob = prob + U_calc_model_order(IMG, temp_SP.N, is_old);
            prob = prob - U_calc_model_order(IMG, temp_SP.N-1, is_old);
        end
    end

    if (prob>max_prob+1e-10 || max_k==-2)
        max_prob = prob;
        max_k = new_k;
    end
end

function top_ok = check_topology(IMG, index, neighborhood)
    [x, y] = get_x_and_y_from_index(index, IMG.xdim);
    %fprintf('x=%d, y=%d, xdim=%d, ydim=%d\n', x, y, IMG.xdim, IMG.ydim);
    top_ok = ( x<=1 || topology_number_label(IMG, x, y, IMG.label(x-1, y), neighborhood)==1 ) && ...
            ( y<=1 || topology_number_label(IMG, x, y, IMG.label(x, y-1), neighborhood)==1 ) && ...
            ( x>=IMG.xdim || topology_number_label(IMG, x, y, IMG.label(x+1, y), neighborhood)==1 ) && ...
            ( y>=IMG.ydim || topology_number_label(IMG, x, y, IMG.label(x, y+1), neighborhood)==1 );
end

function value = topology_number_label(IMG, x, y, check_label, neighborhood)
    neighborhood(2) = y>1 && IMG.label(x, y-1)==check_label;
    neighborhood(4) = x<IMG.xdim && IMG.label(x+1, y)==check_label;
    neighborhood(6) = y<IMG.ydim && IMG.label(x, y+1)==check_label;
    neighborhood(8) = x>1 && IMG.label(x-1, y)==check_label;
    neighborhood(1) = (neighborhood(2) || neighborhood(8)) && (x>1 && y>1 && IMG.label(x-1, y-1)==check_label);
    neighborhood(3) = (neighborhood(2) || neighborhood(4)) && (x<IMG.xdim && y>1 && IMG.label(x+1, y-1)==check_label);
    neighborhood(5) = (neighborhood(4) || neighborhood(6)) && (x<IMG.xdim && y<IMG.ydim && IMG.label(x+1, y+1)==check_label);
    neighborhood(7) = (neighborhood(6) || neighborhood(8)) && (x>1 && y<IMG.ydim && IMG.label(x-1, y+1)==check_label);
    
    index = 0;
    for i=8:-1:1
        index = index*2;
        if neighborhood(i)
            index = index+1;
        end
    end
    value = IMG.T4Table(index+1);
end