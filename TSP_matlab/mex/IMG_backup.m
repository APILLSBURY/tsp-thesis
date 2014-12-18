% =============================================================================
% == IMG.h
% == --------------------------------------------------------------------------
% == An image class to be used with temporal superpixels.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written by Jason Chang and Donglai Wei 06-20-2013
% =============================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Utility Funcitons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IMG = IMG_U_initialize(IMG)
    IMG.new_pos = new_NormalD(2, IMG.hyper.p_theta, IMG.hyper.p_Delta, true);
    IMG.new_app = new_NormalD(IMG.hyper.a_theta, IMG.hyper.a_Delta, true);
    for i=1:numel(IMG.SP)
        if isempty(IMG.SP(i).old) || ~IMG.SP(i).old
            id = IMG.SP(i).UID;
            if (id==0)
                IMG.SP(i) = new_SP(new_pos, new_app, IMG.max_UID, false, IMG.N);
                max_UID = max_UID + 1;
            else
                IMG.SP(i) = new_SP(new_pos, new_app, id, false, IMG.N);
            end
        else
            tmp_pos = new_NormalD(2, IMG.SP(i).v, IMG.hyper.p_theta, IMG.hyper.p_Delta, false);
            tmp_app = new_NormalD(3, IMG.hyper.a_theta, IMG.hyper.a_Delta, false);
            IMG.SP(i) = new_SP(tmp_pos, tmp_app, IMG.SP(i).UID, true, IMG.SP(i).prev_v, IMG.N);
            IMG.SP(i).pos.mean = IMG.SP(i).p_mu;
            IMG.SP(i).app.mean = IMG.SP(i).a_mu;
        end
    end
    for i=numel(IMG.SP)+1:IMG.K
        IMG.SP(i) = new_SP(new_pos, new_app, IMG.max_UID, false, IMG.N);
        IMG.max_UID = IMG.max_UID+1;
    end
    for i=IMG.K+1:IMG.N
        IMG.SP(i) = [];
    end

    % populate the linked lists and pointers
    for i=1:IMG.N
        if (IMG.label(i)>=1)
            IMG.SP(IMG.label(i)) = SP_add_pixel_init(IMG, i, U_check_border_pix(IMG, i));
            IMG.SP(IMG.label(i)) = SP_update_neighbors_add_self(IMG, i);
        end
    end

    for k=1:IMG.K
        IMG.SP(k) = SP_calculate_log_probs(IMG.SP(k));
    end
end


function SP_new = create_SP_new(IMG)
    SP_new = new_SP(IMG.new_pos, IMG.new_app, 0, false, IMG.N);
end


function E = U_calc_energy(IMG)
    E = 0;
    numEmpty = 0;
    for k=1:numel(IMG.SP)
        %if the UID is empty, the SP is empty
        if(~isempty(IMG.SP(k).UID))
            IMG.SP(k) = SP_calculate_log_probs(IMG.SP(k));
            E = E + IMG.SP(k).log_likelihood + U_calc_model_order(IMG, IMG.SP(k).N, IMG.SP(k).old);
        else
            numEmpty = numEmpty + 1;
        end
    end
    E = numEmpty*(create_SP_new(IMG).log_likelihood + U_calc_model_order(IMG, 0, false));
end

function model_order = U_calc_model_order(IMG, size, is_old)
    if size==0
        model_order = 0;
    elseif is_old
        model_order = IMG.log_beta - 0.5*1.837877066409 - 0.5*IMG.log_area_var - 0.5*(IMG.area-size)^2/IMG.area_var;
    else
        model_order = IMG.log_alpha - 0.5*1.837877066409 - 0.5*IMG.log_area_var - 0.5*(IMG.area-size)^2/IMG.area_var;
    end
end

function [x, y] = getXandYfromIndex(index, xdim)
    x = mod(index,(xdim+1))+1;
    y = floor(index/xdim)+1;
end

function border = U_check_border_pix(IMG, index, cur_label)
    [x, y] = getXandYfromIndex(index, IMG.xdim);
    x1 = max(1, x-1);
    y1 = max(1, y-1);
    x2 = min(xdim, x+1);
    y2 = min(ydim, y+1);
    
    %check to see if a label was specified
    if nargin < 3
        cur_label = IMG.label(x, y);
    end

    %check to see if any of the adjacent pixels are different than
    %the current one
    border = (cur_label ~= IMG.label(x1, y) || cur_label ~= IMG.label(x, y1) || ...
        cur_label ~= IMG.label(x2, y) || cur_label ~= IMG.label(x, y2));
end

% --------------------------------------------------------------------------
% -- U_update_border_changed
% --   Updates all border linked lists and pointers for the neighbors of a
% -- changed pixel at index.
% --
% --   parameters:
% --     - index : the index of the recently changed pixel
% --------------------------------------------------------------------------
function IMG = U_update_border_changed(IMG, index)
    [x, y] = getXandYfromIndex(index, IMG.xdim);
    if (x>1)
        IMG = U_update_border_changed_pixel(IMG, index-1);
    end
    if (y>1)
        IMG = U_update_border_changed_pixel(IMG, index-IMG.xdim);
    end
    if (x<IMG.xdim)
        IMG = U_update_border_changed_pixel(IMG, index+1);
    end
    if (y<IMG.ydim)
        IMG = U_update_border_changed_pixel(IMG, index+IMG.xdim);
    end
end

% --------------------------------------------------------------------------
% -- U_update_border_changed_pixel
% --   Updates all border linked lists and pointers at index
% --
% --   parameters:
% --     - index : the index of the recently changed pixel
% --------------------------------------------------------------------------
function IMG = U_update_border_changed_pixel(IMG, index)
    [x, y] = getXandYfromIndex(index, IMG.xdim);
    if (IMG.label(x, y)>=0)
        if (IMG.borders(index) && ~U_check_border_pix(IMG, index))
            IMG.SP(IMG.label(x, y)).borders(index) = false;
            IMG.borders(index) = false;
        elseif (~IMG.borders(index) && U_check_border_pix(IMG, index))
            IMG.borders(index) = true;
            IMG.SP(IMG.label(x, y)).borders(index) = true;
        end
    end
end

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
    [x, y] = getXandYfromIndex(index, IMG.xdim);
    k = IMG.label(x, y);
    if (k>0 && ~neighbors(k))
        neighbors(k) = true;
    end
end

function IMG = U_relabel_SP(IMG, final)
    last_k = 1;
    for k=1:numel(IMG.SP)
        if ~final || IMG.SP(k).N>0
            if (k~=last_k)
                % move the super pixel to the empty one
                IMG.SP(last_k) = IMG.SP(k);
                IMG.SP(k) = [];

                % relabel the pixels
                for i=1:length(IMG.SP(last_k).pixels)
                    if IMG.SP(last_k).pixels(i)
                        [x, y] = getXandYfromIndex(index, IMG.xdim);
                        IMG.label(x, y) = last_k;
                    end
                end
            end
            last_k = last_k + 1;
        end
    end
    IMG.K = last_k-1;
end


function IMG = U_fix_neighbors_self(IMG, k)
    IMG.SP(k).neighbors = [];
    for i=1:length(IMG.SP(k).borders)
        if IMG.SP(k).borders(i)
            IMG.SP(k) = SP_update_neighbors_add_self(IMG.SP(k), IMG.label, i);
        end
    end
end


function IMG = U_fix_neighbors_neighbors(IMG, k, kignore)
    if narg_in < 3
        kignore = 0;
    end
    for i=1:lenth(IMG.SP(k).neighbors)
        k2 = IMG.SP(k).neighbors(i);
        if k2 > 0
            if k2 ~= kignore
                IMG = U_fix_neighbors_self(IMG, k2);
            end
        end
    end
end

% --------------------------------------------------------------------------
% -- U_update_neighbors_rem
% --   Updates all neighbor lists when a pixel is "removed". A subsequent
% -- call to update_neighbors_add should be completed right after this one.
% -- The neighboring label should be changed before calling this function.
% --
% --   parameters:
% --     - old_label : the label of the pixel before it was changed
% --     - index : the index bordering the removed pixel
% --------------------------------------------------------------------------
function IMG = U_update_neighbors_rem(IMG, old_label, index)
    if (old_label>0)
        [x, y] = getXandYfromIndex(index, IMG.xdim);
        if x>1
            llabel = IMG.label(x-1, y);
        else
            llabel = old_label;
        end
        if x<IMG.xdim
            rlabel = IMG.label(x+1, y);
        else
            rlabel = old_label;
        end
        if y>1
            ulabel = IMG.label(x, y-1);
        else
            ulabel = old_label;
        end
        if y<IMG.ydim
            dlabel = IMG.label(x, y+1);
        else
            dlabel = old_label;
        end

        if (llabel~=old_label)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), llabel);
        end
        if (rlabel~=old_label && rlabel ~=llabel)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), rlabel);
        end
        if (ulabel~=old_label && ulabel~=llabel && ulabel~=rlabel)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), ulabel);
        end
        if (dlabel~=old_label && dlabel~=llabel && dlabel~=rlabel && dlabel~=ulabel)
            IMG.SP(old_label) = SP_update_neighbors_label_rem(IMG.SP(old_label), dlabel);
        end

        % update the neighbors' neighbors list. (not a typo!)
        if (llabel~=old_label)
            IMG.SP(llabel) = SP_update_neighbors_label_rem_check(IMG.SP(llabel), IMG.label, x-1, y, old_label);
        end
        if (rlabel~=old_label)
            IMG.SP(rlabel) = SP_update_neighbors_label_rem_check(IMG.SP(rlabel), IMG.label, x+1, y, old_label);
        end
        if (ulabel~=old_label)
            IMG.SP(ulabel) = SP_update_neighbors_label_rem_check(IMG.SP(ulabel), IMG.label, x, y-1, old_label);
        end
        if (dlabel~=old_label)
            IMG.SP(dlabel) = SP_update_neighbors_label_rem_check(IMG.SP(dlabel), IMG.label, x, y+1, old_label);
        end
    end
end

% --------------------------------------------------------------------------
% -- U_update_neighbors_add
% --   Updates all neighbor lists when a pixel is "added". A previous
% -- call to update_neighbors_rem should be completed right before this one.
% -- The neighboring label should be changed before calling this function.
% --
% --   parameters:
% --     - index : the index bordering the removed pixel
% --------------------------------------------------------------------------
function IMG = U_update_neighbors_add(IMG, index)
    [x, y] = getXandYfromIndex(index, IMG.xdim);
    cur_label = IMG.label(x, y);
    if (cur_label>0)
        [x, y] = getXandYfromIndex(index, IMG.xdim);
        if x>1
            llabel = IMG.label(x-1, y);
        else
            llabel = cur_label;
        end
        if x<IMG.xdim
            rlabel = IMG.label(x+1, y);
        else
            rlabel = cur_label;
        end
        if y>1
            ulabel = IMG.label(x, y-1);
        else
            ulabel = cur_label;
        end
        if y<IMG.ydim
            dlabel = IMG.label(x, y+1);
        else
            dlabel = cur_label;
        end

        if (llabel~=cur_label)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), llabel);
        end
        if (rlabel~=cur_label && rlabel ~=llabel)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), rlabel);
        end
        if (ulabel~=cur_label && ulabel~=llabel && ulabel~=rlabel)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), ulabel);
        end
        if (dlabel~=cur_label && dlabel~=llabel && dlabel~=rlabel && dlabel~=ulabel)
            IMG.SP(cur_label) = SP_update_neighbors_label_add(IMG.SP(cur_label), dlabel);
        end

        % update the neighbors' neighbors list. (not a typo!)
        if (llabel~=cur_label)
            IMG.SP(llabel) = SP_update_neighbors_label_add_check(IMG.SP(llabel), IMG.label, x-1, y, cur_label);
        end
        if (rlabel~=cur_label)
            IMG.SP(rlabel) = SP_update_neighbors_label_add_check(IMG.SP(rlabel), IMG.label, x+1, y, cur_label);
        end
        if (ulabel~=cur_label)
            IMG.SP(ulabel) = SP_update_neighbors_label_add_check(IMG.SP(ulabel), IMG.label, x, y-1, cur_label);
        end
        if (dlabel~=cur_label)
            IMG.SP(dlabel) = SP_update_neighbors_label_add_check(IMG.SP(dlabel), IMG.label, x, y+1, cur_label);
        end
    end
end

% --------------------------------------------------------------------------
% -- U_update_neighbors_merge
% --   Updates the neighbor lists for merging two super pixels.
% -- All labels and SP_arr things should be updated *before* calling this.
% --
% --   parameters:
% --     - index : the index bordering the removed pixel
% --------------------------------------------------------------------------
function IMG = U_update_neighbors_merge(IMG, new_label, old_label)
    IMG.SP(new_label).neighbors = IMG.SP(new_label).neighbors + IMG.SP(old_label);
    %The new label can't have itself as its neighbor
    IMG.SP(new_label).neighbors(new_label) = 0;

    % now update all neighboring neighbor lists
    for neighbor_k=1:length(IMG.SP(new_label).neighbors);
        if neighbor_k > 0
            IMG.SP(neighbor_k) = SP_update_neighbors_self(IMG.SP(neighbor_k), IMG.label);
        end
    end
end

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
            SP_update_neighbors_self(IMG.SP(neighbor_k), IMG.label);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Moving Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IMG = Flow_QP2(IMG)
    if (IMG.prev_label~=0)
        %U_find_prev_border();
        %given label z, globally update flow variable of SP using GP approximation of Optical Flow like L2 penalty
        %global variable

        % get the alive stuff
        num_alive = 1;
        alive2all = zeros(IMG.prev_K, 1);
        for k=1:IMG.prev_K
            if (~isempty(IMG.SP(k).N))
                alive2all(num_alive) = k;
                num_alive = num_alive + 1;
            end
        end

        % for computational efficiency
        % 1. Covariance matrix
        if (IMG.alive_dead_changed)
            Syy = zeros(num_alive, num_alive);
            Sxy = zeros(IMG.prev_K, num_alive);

            %obs_uv(1) = obs_u, obs_uv(2) = obs_v
            obs_uv = zeros(num_alive,2);

            IMG.alive_dead_changed = false;
        end

        % populate the observation
        for k1=1:num_alive
            all_k1 = alive2all(k1);
            obs_uv(k1,:) = NormalD_get_mean(IMG.SP(all_k1).pos) - IMG.prev_pos_mean(all_k1) - IMG.SP(all_k1).prev_v;
        end

        % uv_gsl(1) = ugsl, uv_gsl(2) = vgsl
        uv_gsl = Sxy * Syy * obs_uv;

        % copy over the new parameters
        for i=1:IMG.prev_k
            if (~isempty(IMG.SP(i).N))
                % ith element of ugsl
                IMG.SP(i).pos.offset = uv_gsl(i,:) + IMG.SP(i).prev_v;
            end
        end
        IMG.SP_changed(1:IMG.prev_k) = true;
    end
end


% --------------------------------------------------------------------------
% -- link_cost
% --   calculates the energy of linking IMG.SP(ko) (old) to IMG.SP(kn) (new)
% --
% --   parameters:
% --     - ko : the old k
% --     - kn : the new k
% --------------------------------------------------------------------------
function logprob = link_cost(IMG, ko, kn)
    if (ko<IMG.prev_K && kn<IMG.K) % matched SPs
        logprob = SP_log_likelihood_switch_app_prior(IMG.SP(kn), IMG.SP(ko));
        logprob = logprob + U_calc_model_order(IMG, IMG.SP(kn).N, true);
    elseif (ko>=IMG.prev_K) % new SP
        logprob = SP_log_likelihood_switch_app_prior(IMG.SP(kn), create_SP_new(IMG));
        logprob = logprob + U_calc_model_order(IMG, IMG.SP(kn).N, false);
    elseif (kn>=IMG.K) % dead SP
        logprob = SP_log_likelihood_switch_app_prior(create_SP_new(IMG), IMG.SP(ko));
    else
        disp('Cant link a dead SP to a new SP');
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
    [x, y] = getXandYfromIndex(index, IMG.xdim);

    if (new_k<1)
        if IMG.boundary_mask(x,y)
            prob = -inf;
        else
            prob = dummy_log_prob;
        end
    else
        if (new_k > numel(IMG.SP) || isempty(IMG.SP(new_k).N)) % new super pixel
            temp_SP = SP_new;
            is_old = false;
        else
            temp_SP = IMG.SP(new_k);
            is_old = IMG.SP_old(new_k);
        end

        if (add)
            prob = prob + SP_log_likelihood_test_point(temp_SP, IMG.data(index,:), boundary_mask(x, y));
            prob = prob - temp_SP.log_likelihood;

            prob = prob + U_calc_model_order(IMG, temp_SP.N+1, is_old);
            prob = prob - U_calc_model_order(IMG, temp_SP.N, is_old);
        else
            prob = prob + temp_SP.log_likelihood;
            prob = prob - SP_log_likelihood_test_point_rem(temp_SP, IMG.data(index,:), boundary_mask(x, y));

            prob = prob + U_calc_model_order(IMG, temp_SP.N, is_old);
            prob = prob - U_calc_model_order(IMG, temp_SP.N-1, is_old);
        end
    end

    if (prob>max_prob+1e-10 || max_k==-2)
        max_prob = prob;
        max_k = new_k;
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
% --     - nindex : the neighbor to add this point to
% --     - max_prob : the maximum probability of all neighbors
% --     - max_k : the super pixel index of the maximum probability
% --------------------------------------------------------------------------
function [max_prob, max_k] = move_local_calc_delta_MM(IMG, index, new_k, max_prob, max_k)
    prob = SP_log_likelihood_test_point_MM(IMG.SP(new_k), IMG.data(index, :));

    if (prob>max_prob+1e-10 || max_k<0)
        max_prob = prob;
        max_k = new_k;
    end
end


function logprob = move_switch_calc_delta(IMG, oldSP, newSP)
    logprob = SP_log_likelihood_switch_prior(oldSP, newSP);
    logprob = logprob - oldSP.log_likelihood();
    logprob = logprob + U_calc_model_order(IMG, oldSP.N, oldSP.is_old);
    logprob = logprob - U_calc_model_order(IMG, oldSP.N, newSP.is_old);
end



function move_switch_IMG(IMG)
    empty_SPs = false(length(IMG.SP), 1);
    for k=1:length(IMG.SP)
        if (isempty(IMG.SP(k).N) && IMG.SP_old(k))
            empty_SPs(k) = true;
        end
    end

    for k=1:length(IMG.SP)
        if (~isempty(IMG.SP(k).N))
            % if old, check to see a new or unused old
            % if new, check to see if any unused old
            best_k = -1;
            best_energy = 0;
            if (IMG.SP_old(k))
                delta = move_switch_calc_delta(IMG.SP(k), create_SP_new(IMG));
                if (delta > best_energy)
                    best_k = IMG.K;
                    best_energy = delta;
                end
            end

            for test_k=1:length(empty_SPs)
                if empty_SPs(test_k)
                    delta = move_switch_calc_delta(IMG.SP(k), IMG.SP(test_k));
                end
                if (delta > best_energy)
                    best_k = test_k;
                    best_energy = delta;
                end
            end

            % switch with best label
            if (best_k>0 && best_energy>0)
                % change the labels
                for index = 1:length(IMG.SP(k).pixels)
                    if IMG.SP(k).pixels(index)
                        [x, y] = getXandYfromIndex(index, IMG.xdim);
                        IMG.label(x, y) = best_k;
                    end
                end
                IMG.SP_changed(k) = true;
                IMG.SP_changed(best_k) = true;

                % make room for the new one
                if (best_k==IMG.K)
                    if (IMG.K>=IMG.N)
                        disp('Ran out of space!');
                    end
                    if (isempty(IMG.SP(k).N))
                        IMG.SP(k) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, false, IMG.N);
                        IMG.max_UID = IMG.max_UID + 1;
                    end
                    IMG.K = IMG.K + 1;
                else
                    IMG.alive_dead_changed = true;
                end
                SP_merge_with(IMG.SP(best_k), IMG.SP(k), IMG.label);
                
                % delete it if it was a new one
                if (~IMG.SP_old(k))
                    IMG.SP(k) = [];
                else
                    % add it to the search list
                    empty_SPs(k) = true;
                    IMG.alive_dead_changed = true;
                end

                % remove it from the search list
                if (IMG.SP_old(best_k))
                    empty_SPs(best_k) = false;
                end

                % now update the neighbors list
                IMG = U_fix_neighbors_self(IMG, best_k);
                % update the neighbors' neighbors
                IMG = U_fix_neighbors_neighbors(IMG, best_k, IMG.N+1);
            end
        end
    end
end




% --------------------------------------------------------------------------
% -- move_local
% --   finds the optimal local joint move of labels and parameters. Chooses
% -- a super pixel at random and loops through and updates its borders.
% --------------------------------------------------------------------------
function IMG_changed = move_local_IMG(IMG)

    if (IMG.K>IMG.N)
        disp('Ran out of space!');
    end

    changed = false;

    % temporary neighborhood
    neighborhood = false(9,1);

    Nsp = numel(IMG.SP); % number of superpixels

    % choose a random order of super pixels
    perm = randperm(Nsp);

    % MM STEP
    for ki=1:Nsp
        k = perm(ki);

        % find a nonempty super pixel
        if (~isempty(IMG.SP(k).N) && IMG.SP_changed(k))
            IMG.SP_changed(k) = false;
            % loop through borders
            border_length = length(IMG.SP(k).borders);
            i = 1;
            bfs_length = 5 * sum(IMG.SP(k).borders);
            
            %2 wide for x and y values
            bfs_queue = zeros(bfs_length);
            bfs_index = 1;
            
            while (i<=border_length+bfs_length)
                if (i<=border_length)
                    while ~IMG.SP(k).borders(i)
                        i = i+1;
                    end
                    index = i;
                    [x, y] = getXandYfromIndex(index, IMG.xdim);
                else
                    if bfs_queue(i - border_length) == 0
                        break;
                    else
                        index = bfs_queue(i - border_length);
                    end
                    [x, y] = getXandYfromIndex(index, IMG.xdim);
                    k = IMG.label(x, y);
                    if (k>0)
                        i = i+1;
                        continue;
                    end
                end
                i = i+1;

                % check the topology for this pixel
                if (~check_topology(IMG, index, neighborhood))
                    if (k>0)
                        IMG.SP(k).borders(index) = true;
                    end
                else
                    % temporarily remove this data point from the SP
                    max_prob = -inf;
                    max_k = -2;

                    % the current k has to be a possible choice
                    [max_prob, max_k] = move_local_calc_delta(index, IMG.label(x, y), false, max_prob, max_k);

                    if (~boundary_mask(x, y))
                        [max_prob, max_k] = move_local_calc_delta(index, -1, true, max_prob, max_k);
                    end

                    % find which k's we can move to
                    if (x>1 && IMG.label(x-1, y)~=k)
                        [max_prob, max_k] = move_local_calc_delta(index, IMG.label(x-1, y), true, max_prob, max_k);
                    end
                    if (y>1 && IMG.label(x, y-1)~=k)
                        [max_prob, max_k] = move_local_calc_delta(index, IMG.label(x, y-1), true, max_prob, max_k);
                    end
                    if (x<xdim && IMG.label(x+1, y)~=k)
                        [max_prob, max_k] = move_local_calc_delta(index, IMG.label(x+1, y), true, max_prob, max_k);
                    end
                    if (y<ydim && IMG.label(x, y+1)~=k)
                        [max_prob, max_k] = move_local_calc_delta(index, IMG.label(x, y+1), true, max_prob, max_k);
                    end

                    if (max_k~=k)
                        changed = true;
                        % update the labels... it moves from k->max_k
                        IMG.label(x, y) = max_k;
                        if (max_k>=IMG.K)
                            % creating a new one
                            if (IMG.K>=N)
                                disp('Ran out of space!');
                            end
                            IMG.SP(k) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, false, IMG.N);
                            IMG.max_UID = IMG.max_UID + 1;
                            IMG.K = IMG.K + 1;
                        end

                        % update the neighbors lists
                        IMG = U_update_neighbors_rem(IMG, k, index);
                        IMG = U_update_neighbors_add(IMG, index);

                        % set the correct SP_changed variables of all neighbors
                        if (k<1)
                            IMG.SP_changed(max_k) = true;
                        else
                            for neighbor_k=1:numel(IMG.SP(k).neighbors)
                                if IMG.SP(k).neighbors(neighbor_k) > 0
                                    IMG.SP_changed(neighbor_k) = true;
                                end
                            end
                            if (max_k>0)
                                for neighbor_k=1:numel(IMG.SP(max_k).neighbors)
                                    if IMG.SP(max_k).neighbors(neighbor_k) > 0
                                        IMG.SP_changed(neighbor_k) = true;
                                    end
                                end
                            end
                        end


                        % update all border lists for neighbors
                        IMG = U_update_border_changed(IMG, index);
                        if (k>0)
                            SP_rem_pixel(IMG.SP(k), IMG.data(index, :), index);
                        end

                        if (k>0 && isempty(IMG.SP(k).N))
                            if (IMG.SP_old(k))
                                IMG.alive_dead_changed = true;
                            else
                                IMG.SP(k) = [];
                                IMG.SP_changed(k) = false;
                            end
                        end

                        % add this point to the maximum SP
                        if (max_k>=0)
                            SP_add_pixel(IMG.SP(max_k), IMG.data(index, :), U_check_border_pix(IMG, index), boundary_mask(x, y));
                        end
                    elseif (k>0)
                        IMG.SP(k).borders(index) = true;
                    end
                    if (x>1 && IMG.label(x-1, y)<1 && bfs_index <= bfs_length)
                        bfs_queue(bfs_index) = index-1;
                        bfs_index = bfs_index + 1;
                    end
                    if (y>1 && IMG.label(x, y-1)<1 && bfs_index <= bfs_length)
                        bfs_queue(bfs_index) = index-xdim;
                        bfs_index = bfs_index + 1;
                    end
                    if (x<xdim && IMG.label(x+1, y)<1 && bfs_index <= bfs_length)
                        bfs_queue(bfs_index) = index+1;
                        bfs_index = bfs_index + 1;
                    end
                    if (y<ydim && IMG.label(x, y+1)<1 && bfs_index <= bfs_length)
                        bfs_queue(bfs_index) = index+xdim;
                        bfs_index = bfs_index + 1;
                    end
                end
            end
        end
    end
    if (~changed)
        IMG.SP_changed(1:IMG.K) = false;
    end
    IMG_changed.IMG = IMG;
    IMG_changed.changed = changed;
end


% --------------------------------------------------------------------------
% -- move_merge_calc_delta
% --   calculates the probability of assigning the pixel at index to the
% -- cluster of nindex (neighbor index).
% --
% --   parameters:
% --     - index : the new point to add
% --     - nindex : the neighbor to add this point to
% --------------------------------------------------------------------------
function prob = move_merge_calc_delta(IMG, k, merge_k)
    prob = SP_log_likelihood_test_merge(IMG.SP(k), IMG.SP(merge_k)) + IMG.SP(merge_k).log_likelihood_empty;
    prob = prob - IMG.SP(k).log_likelihood + IMG.SP(merge_k).log_likelihood;

    prob = prob + U_calc_model_order(IMG, IMG.SP(k).N + IMG.SP(merge_k).N, IMG.SP_old(k));
    prob = prob + U_calc_model_order(IMG, 0, IMG.SP_old(merge_k));
    prob = prob - U_calc_model_order(IMG, IMG.SP(k).N, IMG.SP_old(k));
    prob = prob - U_calc_model_order(IMG, IMG.SP(merge_k).N, IMG.SP_old(merge_k));
end

% --------------------------------------------------------------------------
% -- move_split_calc_delta
% --   calculates the change in energy for
% -- (k1 U new_k1) && (k2 U new_k2) - (k1 U new_k1 U new_k) && (k2)
% --
% --   parameters:
% --     - SP1 : the SP that originates the splitting
% --     - SP2 : the SP to split to
% --     - new_SP1 : temporary SP that contains pixels that will go in k1
% --     - new_SP2 : temporary SP that contains pixels that will go in k2
% --     - SP1_old : indicates if SP1 is an old SP
% --     - SP2_old : indicates if SP2 is an old SP
% --------------------------------------------------------------------------
function prob = move_split_calc_delta(IMG, SP1, SP2, new_SP1, new_SP2, SP1_old, SP2_old)
    prob = SP_log_likelihood_test_merge(SP1, new_SP1);
    prob = prob + SP_log_likelihood_test_merge(SP2, new_SP2);
    prob = prob - SP_log_likelihood_test_merge(SP1, new_SP1, new_SP2);
    prob = prob - SP2.log_likelihood;

    % split
    prob = prob + U_calc_model_order(IMG, SP1.N + new_SP1.N, SP1_old);
    prob = prob + U_calc_model_order(IMG, SP2.N + new_SP2.N, SP2_old);

    % not split
    prob = prob - U_calc_model_order(IMG, SP1.N + new_SP1.N + new_SP2.N, SP1_old);
    prob = prob - U_calc_model_order(IMG, SP2.N, SP2_old);
end

% --------------------------------------------------------------------------
% -- move_merge
% --   finds the optimal merge between all pairs of super pixels
% --------------------------------------------------------------------------
function IMG = move_merge_IMG(IMG)
    % choose a random order of super pixels
    Nsp = numel(IMG.SP);
    perm = randperm(Nsp);

    neighbors = false(Nsp, 1);

    for ki=1:Nsp
        k = perm(ki);

        % find a nonempty super pixel
        if ~isempty(IMG.SP(k).N)
            % find all bordering super pixels
            
            neighbors = U_find_border_SP(IMG, k, neighbors);

            max_E = -inf;
            max_k = -1;

            % loop through all neighbors
            for merge_k=1:length(neighbors)
                if neighbors(merge_k)
                    new_E = move_merge_calc_delta(IMG, k, merge_k);
                    if (new_E > max_E || max_k==-1)
                        max_E = new_E;
                        max_k = merge_k;
                    end

                    % fix neighbors for the next super pixel
                    neighbors(merge_k) = false;
                end
            end

            % merge if it increases energy
            if (max_E>0)
                % change the labels
                for index=1:length(IMG.SP(max_k).pixels)
                    if IMG.SP(max_k).pixels(index)
                        [x, y] = getXandYfromIndex(index, IMG.xdim);
                        IMG.label(x, y) = k;
                    end
                end
                IMG.SP_changed(k) = true;
                IMG.SP_changed(max_k) = true;

                SP_merge_with(IMG.SP(k), IMG.SP(max_k), IMG.label);
                if (~IMG.SP_old(max_k))
                    IMG.SP(max_k) = [];
                else
                    IMG.alive_dead_changed = true;
                end
            end
        end
    end
end


function IMG = move_IMG(IMG)
   % choose a random order of super pixels
   Nsp = numel(IMG.SP);
   perm = randperm(Nsp);

   pre_K = Nsp;
   %split_thres = floor(IMG.N/IMG.K);

   energies = zeros(Nsp, 1);
   max_energy = -inf;
   min_energy = inf;
   mean_area = 0;
   
   for k=1:pre_K
      temp_energy = (IMG.SP(k).log_likelihood + U_calc_model_order(IMG, IMG.SP(k).N, IMG.SP_old(k))) / IMG.SP(k).N;
      energies(k) = temp_energy;
      max_energy = max(max_energy, temp_energy);
      min_energy = min(min_energy, temp_energy);
      mean_area = mean_area + IMG.SP(k).N;
   end
   mean_area = mean_area / pre_K;
   threshold = min_energy + (max_energy-min_energy)*0.2;

   for ki=1:pre_K
      k = perm(ki);
      if (sum(IMG.SP(k).neighbors) > 0 && (IMG.SP(k).N>mean_area || energies(k) < threshold) && IMG.SP_changed(k))
         move_split_SP(k);
      end
   end
end



function IMG = move_local_SP(IMG, k)
    if k > 0 && k <= numel(IMG.SP) && ~isempty(IMG.SP(k).N)
        % loop through borders
        int length = IMG.SP(k)->borders.getLength();
        for index=1:length(IMG.SP(k).borders)
            if IMG.SP(k).borders(index)
                if isempty(IMG.SP(k).N)
                    break;
                end

                max_prob = -10^10;
                max_k = -1;
                [x, y] = getXandYfromIndex(index, IMG.xdim);

                % the current k has to be a possible choice
                [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x, y), false, max_prob, max_k);
                % a new k is also a possible choice
                [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.K, true, max_prob, max_k);

                % find which k's we can move to
                if (x>1 && IMG.label(x-1, y)~=k)
                    [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x-1, y), true, max_prob, max_k);
                end
                if (y>1 && IMG.label(x, y-1)~=k)
                    [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x, y-1), true, max_prob, max_k);
                end
                if (x<xdim && IMG.label(x+1, y)~=k)
                    [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x+1, y), true, max_prob, max_k);
                end
                if (y<ydim && IMG.label(x, y+1)~=k)
                    [max_prob, max_k] = move_local_calc_delta(IMG, index, IMG.label(x, y+1), true, max_prob, max_k);
                end


                if (max_k~=k)
                    % update the labels
                    IMG.label(x, y) = max_k;
                    if (max_k>=IMG.K)
                        IMG.SP(k) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, false, IMG.N);
                        IMG.max_UID = IMG.max_UID + 1;
                        IMG.K = IMG.K + 1;
                    end

                    % update all border lists for neighbors
                    IMG = U_update_border_changed(IMG, index);

                    SP_rem_pixel(IMG.SP(k), IMG.data(index,:), IMG.boundary_mask(x, y));
                    if isempty(IMG.SP(k).N)
                        if (IMG.SP_old(k))
                            IMG.alive_dead_changed = true;
                        else
                            IMG.SP(k) = [];
                        end
                    end

                    % add this point to the maximum SP
                    SP_add_pixel(IMG.SP(max_k), IMG.data(index,:), U_check_border_pix(IMG, index), IMG.boundary_mask(x, y));
                else
                    IMG.SP(k).borders(index) = true;
                end
            end
        end
    end
end


function neighbors = U_find_border_SP(IMG, k, neighbors)
    for index=1:length(IMG.SP(k).borders)
        if IMG.SP(k).borders(index)
            [x, y] = getXandYfromIndex(index, IMG.xdim);

            if (x>1 && IMG.label(x-1, y)~=k)
                neighbors = U_update_neighbor_list(IMG, neighbors, index-1);
            end
            if (y>1 && IMG.label(x, y-1)~=k)
                neighbors = U_update_neighbor_list(IMG, neighbors, index-IMG.xdim);
            end
            if (x<IMG.xdim && IMG.label(x+1, y)~=k)
                neighbors = U_update_neighbor_list(IMG, neighbors, index+1);
            end
            if (y<IMG.ydim && IMG.label(x, y+1)~=k)
                neighbors = U_update_neighbor_list(IMG , neighbors, index+IMG.xdim);
            end
        end
    end
end


function [max_E, max_k] = move_merge_SP_propose(IMG, k, neighbors, max_E, max_k)
    U_find_border_SP(k,neighbors);
    % loop through all neighbors
    for merge_k=1:length(neighbors)
        if neighbors(merge_k)
            % calculate the energy
            tmp_E = move_merge_calc_delta(IMG, k, merge_k);
            if(tmp_E>max_E)
                max_E = tmp_E;
                max_k = merge_k;
            end
            % fix neighbors for the next super pixel
            neighbors(merge_k) = false;
        end
    end
    if (max_k==-1)
        disp('find_bestnb: No neighbour found');
    end
end


void IMG::move_merge_SP_propose_region(int k,arr(bool) neighbors,linkedList<int> &check_labels,double &max_E,int &max_k)
{
   linkedList<int> neighborsLL;
   U_find_border_SP(k,neighbors,neighborsLL);

   % loop through all neighbors
   linkedListNode<int>* node = neighborsLL.getFirst();
   double tmp_E = 0;
   int merge_k;
   while (node != NULL)
   {
      merge_k = node->getData();
      %printf("inside %d,%d,%d\n",merge_k,k,check_labels.getLength());
      if(U_FindArray(merge_k,check_labels)!=-1)
      {
         % calculate the energy
         tmp_E = move_merge_calc_delta(k, merge_k);
         if(tmp_E>max_E)
         {
               max_E = tmp_E;
               max_k = merge_k;
         }
      }
      % fix neighbors for the next super pixel
      neighbors[merge_k] = false;
      node = node->getNext();
   }
   if (max_k==-1)
   {
      neighborsLL.print();
      mexPrintf("k=%d\n", debug);
      mexErrMsgTxt("find_bestnb_region: No neighbour found\n");
   }
}
void IMG::move_merge_SP(int k,arr(bool) neighbors)
{
   % find a nonempty super pixel
   if (IMG.SP(k)!=NULL && !(IMG.SP(k)->isempty()))
   {
      double max_E = -DBL_MAX;
      int max_k = -1;
      move_merge_SP_propose( k, neighbors,max_E,max_k);
      % merge if it increases energy
      if (max_E>0)
      {
         % change the labels
         linkedListNode<int>* node = IMG.SP(max_k]->pixels.getFirst();
         while (node != NULL)
         {
            int index = node->getData();
            IMG.label[index] = k;
            node = node->getNext();
         }
         IMG.SP(k)->merge_with(IMG.SP(max_k], IMG.label, border_ptr, xdim, ydim);
         if (IMG.SP_old[max_k])
         {
            delete IMG.SP(max_k];
            IMG.SP(max_k] = NULL;
         }
      }
   }
}






void IMG::move_split_SP_propose(int index, int num_SP, int option, double &max_E,int &ksplit,int* new_ks)
{
   int num_iter = 5,SP_bbox[4]={xdim,0,ydim,0};
   % 1. Kmeans++
   bool broken = U_Kmeans_plusplus(index, SP_bbox, num_SP, num_iter, option);

   if (broken)
   {
      max_E = -DBL_MAX;
      ksplit = -2;
      return;
   }

   % 2. create two new SP from the old one
   % label matrix is equivalent to SP.pixels
   % we will mess around label matrix for new proposal
   % and recover it from SP.pixels if it doesn't work out
   % merge small connected component in the label matrix
   U_connect_newSP(SP_bbox, num_SP);

   % new super pixels in K and K+1... old super pixel in index
   IMG.SP(index]->empty(false);
   double E;
   if (option==-1)
   {
      % (index, new);
      E = move_split_calc_delta(IMG.SP(index], SP_new, IMG.SP(k), IMG.SP(K+1], IMG.SP_old[index], false);
      if (E>max_E)
      {
         max_E = E;
         new_ks[0] = K;
         new_ks[1] = K+1;
         ksplit = -1;
      }
      % (new, index)
      E = move_split_calc_delta(IMG.SP(index], SP_new, IMG.SP(K+1], IMG.SP(k), false, IMG.SP_old[index]);
      if (E>max_E)
      {
         max_E = E;
         new_ks[0] = K+1;
         new_ks[1] = K;
         ksplit = -1;
      }
      % (index, old_empty) && (old_empty, index)
      for (int ktest=0; ktest<K; ktest++) if (ktest!=index && IMG.SP(ktest]!=NULL && IMG.SP(ktest]->get_N()==0)
      {
         E = move_split_calc_delta(IMG.SP(index], IMG.SP(ktest], IMG.SP(k), IMG.SP(K+1], IMG.SP_old[index], IMG.SP_old[ktest]);
         if (E>max_E)
         {
            max_E = E;
            new_ks[0] = K;
            new_ks[1] = K+1;
            ksplit = ktest;
         }

         E = move_split_calc_delta(IMG.SP(index], IMG.SP(ktest], IMG.SP(K+1], IMG.SP(k), IMG.SP_old[index], IMG.SP_old[ktest]);
         if (E>max_E)
         {
            max_E = E;
            new_ks[0] = K+1;
            new_ks[1] = K;
            ksplit = ktest;
         }
      }

   }
   else
   {
      % for refine_move
      % split into IMG.SP(index] and IMG.SP(option]
      double E_1_K1;
      double E_1_K2;
      if (!IMG.SP_old[index] && !IMG.SP_old[option])
      {
         arr(double) SP1_total_pos = IMG.SP(K-2]->get_total_pos();
         double N1 = IMG.SP(K-2]->get_N();
         arr(double) SP2_total_pos = IMG.SP(K-1]->get_total_pos();
         double N2 = IMG.SP(K-1]->get_N();
         arr(double) SPK1_total_pos = IMG.SP(k)->get_total_pos();
         double NK1 = IMG.SP(k)->get_N();
         arr(double) SPK2_total_pos = IMG.SP(K+1]->get_total_pos();
         double NK2 = IMG.SP(K+1]->get_N();

         E_1_K1 = 0;
         E_1_K2 = 0;
         for (int d=0; d<2; d++)
         {
            double temp = SPK1_total_pos[d]/NK1 - SP1_total_pos[d]/N1;
            E_1_K1 += temp*temp;
            temp = SPK2_total_pos[d]/NK2 - SP2_total_pos[d]/N2;
            E_1_K1 += temp*temp;

            temp = SPK1_total_pos[d]/NK1 - SP2_total_pos[d]/N2;
            E_1_K2 += temp*temp;
            temp = SPK2_total_pos[d]/NK2 - SP1_total_pos[d]/N1;
            E_1_K2 += temp*temp;
         }

         if (E_1_K1 < E_1_K2)
         {
            %if ((index==177 && option==493) || (index==493 && option==177)){mexPrintf("1) E1K1=%e, E1K2=%e\n", E_1_K1, E_1_K2);drawnow();}
            new_ks[0] = K;
            new_ks[1] = K+1;
            max_E = move_split_calc_delta(IMG.SP(index], IMG.SP(option], IMG.SP(k), IMG.SP(K+1], IMG.SP_old[index], IMG.SP_old[option]);
         }
         else
         {
            %if ((index==177 && option==493) || (index==493 && option==177)){mexPrintf("2) E1K1=%e, E1K2=%e\n", E_1_K1, E_1_K2);drawnow();}
            new_ks[0] = K+1;
            new_ks[1] = K;
            max_E = move_split_calc_delta(IMG.SP(index], IMG.SP(option], IMG.SP(K+1], IMG.SP(k), IMG.SP_old[index], IMG.SP_old[option]);
         }
      }
      else
      {
         E_1_K1 = move_split_calc_delta(IMG.SP(index], IMG.SP(option], IMG.SP(k), IMG.SP(K+1], IMG.SP_old[index], IMG.SP_old[option]);
         E_1_K2 = move_split_calc_delta(IMG.SP(index], IMG.SP(option], IMG.SP(K+1], IMG.SP(k), IMG.SP_old[index], IMG.SP_old[option]);

         if (E_1_K1 > E_1_K2)
         {
            new_ks[0] = K;
            new_ks[1] = K+1;
            max_E = E_1_K1;
         }
         else
         {
            new_ks[0] = K+1;
            new_ks[1] = K;
            max_E = E_1_K2;
         }
      }
   }
}



void IMG::move_split_SP(int index)
{
   int num_SP = 2;
   if (IMG.SP(index]!=NULL && !(IMG.SP(index]->isempty()) && IMG.SP(index]->get_N() > num_SP)
   {

      int* new_ks = new int[num_SP];
      memset(new_ks,-1,sizeof(int)*num_SP);
      int ksplit = -1;
      double max_E = -DBL_MAX;

      %only work for num_SP == 2, so far
      move_split_SP_propose(index, num_SP, -1, max_E, ksplit, new_ks);

      % update
      if(max_E>0)
      {
         %mexPrintf("split: %f,%d,%d\n",max_E,index,new_ks[0]);
         % update the labels first
         linkedListNode<int>* node = IMG.SP(new_ks[0]]->pixels.getFirst();
         while (node != NULL)
         {
            IMG.label[node->getData()] = index;
            node = node->getNext();
         }

         % merge the super pixels
         IMG.SP(index]->merge_with(IMG.SP(new_ks[0]], IMG.label, border_ptr, xdim, ydim);
         IMG.SP(new_ks[0]]->empty();
         delete IMG.SP(new_ks[0]];
         IMG.SP(new_ks[0]] = NULL;

         IMG.SP_changed[index] = true;
         if (ksplit<0)
         {
            % splitting into a new super pixel
            IMG.SP_changed[new_ks[1]] = true;
            % move it to the right spot
            if (new_ks[1]!=K)
            {
               node = IMG.SP(new_ks[1]]->pixels.getFirst();
               while (node != NULL)
               {
                  IMG.label[node->getData()] = K;
                  node = node->getNext();
               }
               if (IMG.SP(k)!=NULL)
                  mexPrintf("%d\n",IMG.SP(k)->get_N());
               IMG.SP(k) = IMG.SP(new_ks[1]];
               IMG.SP_old[K] = IMG.SP_old[new_ks[1]];
               IMG.SP(new_ks[1]] = NULL;
            }
            K++;
            max_UID++;
         }
         else
         {
            % splitting into an old super pixel
            IMG.SP_changed[ksplit] = true;
            node = IMG.SP(new_ks[1]]->pixels.getFirst();
            while (node != NULL)
            {
               IMG.label[node->getData()] = ksplit;
               node = node->getNext();
            }
            IMG.SP(ksplit]->merge_with(IMG.SP(new_ks[1]], IMG.label, border_ptr, xdim, ydim);
            IMG.SP(new_ks[1]]->empty();
            delete IMG.SP(new_ks[1]];
            IMG.SP(new_ks[1]] = NULL;
            IMG.alive_dead_changed = true;
         }
      }
      else if (ksplit!=-2)
      {
         %Recover previous IMG.label
         linkedListNode<int> *node = IMG.SP(k)->pixels.getFirst();
         while (node!=NULL)
         {
            IMG.label[node->getData()] = index;
            node = node->getNext();
         }
         node = IMG.SP(K+1]->pixels.getFirst();
         while (node!=NULL)
         {
            IMG.label[node->getData()] = index;
            node = node->getNext();
         }

         IMG.SP_changed[index] = false;

         for( int i = 0; i < num_SP; i++ )
         {
            IMG.SP(index]->merge_with(IMG.SP(K+i], IMG.label, border_ptr, xdim, ydim);
            IMG.SP(K+i]->empty();
            delete IMG.SP(K+i];
            IMG.SP(K+i] = NULL;
         }
      }
      else
      {
         IMG.SP_changed[index] = false;
      }

      delete[] new_ks;
   }
}






% this is split stuff!
double IMG::U_dist(int index1 , int index2)
{
   double dist = 0;
   int ind1=index1*5,ind2=index2*5;
   %printf("dis:%d,%d\n",ind1,ind2);
   %x,y
   %arr(double) pos_Delta = new_pos.get_Delta();
   %arr(double) app_Delta = new_app.get_Delta();
   for( int n = 0; n < 2; n++ )
   {
     dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
     % dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
   }
   %L,a,b
   for( int n = 2; n < 5; n++ )
   {
      dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
      %dist += (data[ind1+n]-data[ind2+n])*(data[ind1+n]-data[ind2+n]);
   }
      %printf("wa: %f\n",(dist/10000.0));
   return dist/10000.0;
}
double IMG::U_dist(int index1 , double* center)
{
   double dist = 0;
   int ind1=index1*5;

   %arr(double) pos_Delta=this->new_pos.get_Delta();
   %arr(double) app_Delta=this->new_app.get_Delta();

   for( int n = 0; n < 2; n++ ){
      dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
      %dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
   }
   %L,a,b
   for( int n = 2; n < 5; n++ ){
      dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
      %dist += (data[ind1+n]-center[n])*(data[ind1+n]-center[n]);
   }
      %printf("wa2 : %f\n",(dist/10000.0));
      return dist/10000.0;
}
bool IMG::U_Kmeans_plusplus(int index, int *bbox, int num_SP, int numiter, int index2)
{
   bool broken = false;

   if (num_SP!=2)
      mexErrMsgTxt("Trying to split into more than 2");
   %void IMG::U_Kmeans_plusplus( int num_SP, int numiter){
   int num_pix = IMG.SP(index]->get_N();
   %int num_pix = width*height;
   double* distvec = new double[num_pix];
   int* klabels = new int[num_pix];
   int* SP_sz = new int[num_SP];
   double * center = new double[num_SP * 5];
   double dist,tmp_dist;
   int change;

   %1. kmeans ++ initialization
   %first cener pt

   int seed;


   linkedListNode<int>* node = IMG.SP(index]->pixels.getFirst();
   int count;
   int tmp_pos;
   while (node!=NULL)
   {
      %distvec[n++] = U_dist(seed,node->getData());
      %distvec[n++] = U_dist(node->getData(), center);
      % get the Bounding Box of IMG.SP(index]
      % xmin,xmax,ymin,ymax
      % x
      tmp_pos = node->getData()%xdim;
      if(tmp_pos <bbox[0])
         bbox[0] = tmp_pos;
      if(tmp_pos >bbox[1])
         bbox[1] = tmp_pos;
      % y
      tmp_pos = node->getData()/xdim;
      if(tmp_pos <bbox[2])
         bbox[2] = tmp_pos;
      if(tmp_pos >bbox[3])
         bbox[3] = tmp_pos;
      node = node->getNext();
   }

   bool old_split1 = IMG.SP_old[index];
   bool old_split2 = (index2!=-1 && IMG.SP_old[index2]);
   if (!old_split1 && !old_split2) % both new
   {
      % first is new
      seed = IMG.SP(index]->pixels.dataAt(rand()%num_pix);
      if(IMG.label[seed]!=index)
      {
         mexPrintf("IMG.SP(index]->get_N()=%d\n", IMG.SP(index]->get_N());
         mexErrMsgTxt("inconsistency about cluster label\n");
      }
      for( int n = 0; n < 5; n++ )
         center[n] = data[seed*5+n];

      node = IMG.SP(index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         distvec[count++] = U_dist(node->getData(), center);
         node = node->getNext();
      }

      % second is new
      seed = IMG.SP(index]->pixels.dataAt(U_randmult(distvec,num_pix));
      if(IMG.label[seed]!=index)
         mexErrMsgTxt("inconsistency about cluster label\n");
      for( int m = 0; m < 5; m++ )
         center[5+m] = data[seed*5+m];
   }
   else if (old_split1 && !old_split2) % old new
   {
      arr(double) mean_pos = IMG.SP(index]->get_mean_pos();
      arr(double) mean_app = IMG.SP(index]->get_mean_app();
      center[0] = mean_pos[0];
      center[1] = mean_pos[1];
      center[2] = mean_app[0];
      center[3] = mean_app[1];
      center[4] = mean_app[2];

      node = IMG.SP(index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         distvec[count++] = U_dist(node->getData(), center);
         node = node->getNext();
      }

      % second is new
      seed = IMG.SP(index]->pixels.dataAt(U_randmult(distvec,num_pix));
      if(IMG.label[seed]!=index)
         mexErrMsgTxt("inconsistency about cluster label\n");
      for( int m = 0; m < 5; m++ )
         center[5+m] = data[seed*5+m];
   }
   else if (!old_split1 && old_split2) % new old
   {
      arr(double) mean_pos = IMG.SP(index2]->get_mean_pos();
      arr(double) mean_app = IMG.SP(index2]->get_mean_app();
      center[5] = mean_pos[0];
      center[6] = mean_pos[1];
      center[7] = mean_app[0];
      center[8] = mean_app[1];
      center[9] = mean_app[2];

      node = IMG.SP(index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         distvec[count++] = U_dist(node->getData(), center+5);
         node = node->getNext();
      }

      % first is new
      seed = IMG.SP(index]->pixels.dataAt(U_randmult(distvec,num_pix));
      if(IMG.label[seed]!=index)
         mexErrMsgTxt("inconsistency about cluster label\n");
      for( int m = 0; m < 5; m++ )
         center[m] = data[seed*5+m];
   }
   else % old old
   {
      arr(double) mean_pos = IMG.SP(index]->get_mean_pos();
      arr(double) mean_app = IMG.SP(index]->get_mean_app();
      center[0] = mean_pos[0];
      center[1] = mean_pos[1];
      center[2] = mean_app[0];
      center[3] = mean_app[1];
      center[4] = mean_app[2];

      mean_pos = IMG.SP(index2]->get_mean_pos();
      mean_app = IMG.SP(index2]->get_mean_app();
      center[5] = mean_pos[0];
      center[6] = mean_pos[1];
      center[7] = mean_app[0];
      center[8] = mean_app[1];
      center[9] = mean_app[2];
   }


   %2. kmeans ++ iterations
   change = 0;
   for( int itr = 0; itr < numiter; itr++ ){
      memset(distvec, 0x7f, sizeof(double)*num_pix);
      node = IMG.SP(index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {

         for( int i = 0; i < num_SP; i++ ){
            tmp_dist = U_dist(node->getData(),&center[i*5]);
            %				printf("%d,%d,%f\n",n,node->getData(),tmp_dist);
            if(tmp_dist < distvec[count] ){
               distvec[count] = tmp_dist;
               klabels[count]  = i;
               change = 1;
            }
         }
         node = node->getNext();
         count++;
      }

      if (change ==0){
         %no change happened... Kmeans totally stuck
         break;
      }
      memset(center,0,sizeof(double)*num_SP*5);
      memset(SP_sz,0,sizeof(int)*num_SP);


      node = IMG.SP(index]->pixels.getFirst();
      count = 0;
      while (node!=NULL)
      {
         %printf("dist %d,%d,%d,%d\n",r,c,counter,klabels[counter]);
         % klabels[n]==0 && old_split1 then don't update
         % klabels[n]==1 && old_split2 then don't update
         % if (klabels[n]!=old_split1-1 && klabels[n]!=old_split2-1)
         if ((klabels[count]!=0 || !old_split1) && (klabels[count]!=1 || !old_split2))
            for( int j = 0; j < 5; j++ )
               center[klabels[count]*5+j] += data[(node->getData())*5+j];
         SP_sz[klabels[count]] += 1;
         count++;
         node = node->getNext();
      }

      double inv;
      for( int k = 0; k < num_SP; k++ )
      {
         if (SP_sz[k]>0)
         {
            if ((k==0 && !old_split1) || (k==1 && !old_split2))
            {
               inv = 1.0/double(SP_sz[k]);
               for( int kk = 0; kk < 5; kk++ )
                  center[k*5+kk] *= inv;
            }
         }
         else
         {
            if (old_split1 || old_split2)
               broken = true;
            else
               mexErrMsgTxt("one cluster removed... shouldn't happen\n");
         }
      }
   }

   if (!broken)
   {
      % change label accordingly
      node = IMG.SP(index]->pixels.getFirst();
      count = 0;
      while (node != NULL)
      {
         IMG.label[node->getData()] = K+klabels[count++];
         node = node->getNext();
      }
   }

   delete[] distvec;
   delete[] SP_sz;
   delete[] center;
   delete[] klabels;
   return broken;
}



% KMerge Version 2: grow region

void IMG::U_connect_newSP(int *bbox,int num_SP){
   int b, c;
   int min_x = bbox[0];
   int min_y = bbox[2];
   int max_x = bbox[1];
   int max_y = bbox[3];
   int tmp_xdim = max_x-min_x+1;
   int tmp_ydim = max_y-min_y+1;
   int tmp_N = tmp_ydim*tmp_xdim;
   %printf("%d,%d,%d,%d\n",bbox[0],bbox[1],bbox[2],bbox[3]);
   const int dx4[4] = {-1,  0,  1,  0};
   const int dy4[4] = { 0, -1,  0,  1};

   arr(int) nlabels = allocate_memory<int>(tmp_N,-1);
   int oindex = 0, adjlabel = 0, label_count =0;

   int base_ind = min_y*xdim+min_x,tmp_ind,tmp_x,tmp_y,new_ind,new_nind,new_x,new_y;
   linkedListNode<int>* node;
   linkedList<int> check_labels;
   for(int j =0;j<num_SP;j++)
      check_labels.addNodeEnd(K+j);

   for( int j = 0; j < tmp_ydim; j++ )%tmp_ydim
   {
      for( int k = 0; k < tmp_xdim; k++ )%tmp_xdim
      {

         tmp_ind = base_ind+k;
         if(U_FindArray(IMG.label[tmp_ind],check_labels)!=-1 && nlabels[oindex]==-1)
         {
         %printf("do it %d,%d,%d,%d,%d,%d\n",j,k,tmp_ind,base_ind,min_y*xdim,min_x);
            % first find an unlabeled pixel in labels
            nlabels[oindex] = label_count;
            %--------------------
            % Start a new segment of index oindex
            %--------------------
            if(IMG.SP(K + label_count]!= NULL)
               mexErrMsgTxt("SP should be null...\n");
            %mexErrMsgTxt(strcat(,"th SP should be null...\n"));
            IMG.SP(K + label_count] = new SP(new_pos, new_app, max_UID);
            %can't decide whether it will be border yet
            IMG.SP(K + label_count]->add_pixel(data, tmp_ind, false, pixel_ptr[tmp_ind], border_ptr[tmp_ind], boundary_mask[tmp_ind]);

            node = IMG.SP(K + label_count]->pixels.getFirst();
            %int cc=1;
            while (node!=NULL)
            {
               tmp_x = node->getData()%xdim;
               tmp_y = node->getData()/xdim;
               for( int n = 0; n < 4; n++ )
               {
                  new_x = tmp_x+dx4[n];
                  new_y = tmp_y+dy4[n];
                  if( (new_x >= min_x && new_x <= max_x) && (new_y >= min_y && new_y <= max_y) )
                  {
                     new_ind = new_y*xdim + new_x;
                     new_nind = (new_x-min_x)+(new_y-min_y)*tmp_xdim;
                     if(nlabels[new_nind] ==-1 && IMG.label[new_ind] == IMG.label[tmp_ind] )
                     {
                        % should always update labels before adding pixel, otherwise
                        % U_check_border_pix will be wrong!
                        nlabels[new_nind] = label_count;
                        IMG.label[new_ind] = K+label_count;
                        IMG.SP(K + label_count]->add_pixel(data, new_ind, false, pixel_ptr[new_ind], border_ptr[new_ind], boundary_mask[new_ind]);
                     }
                  }
               }
               node = node->getNext();
            }
            IMG.label[tmp_ind] = K+label_count;

            label_count++;
         }
         oindex ++;
      }
      base_ind += xdim;
   }

   % Now is the time to clean up the border pixels and neighbor ids
   for (int k=K; k<K+label_count; k++){
      IMG.SP(k)->fix_borders(IMG.label, border_ptr, xdim, ydim);
      U_fix_neighbors_self(k);
   }

   if (label_count>num_SP)
   {
      %printf("oops... kmeans gives %d connected components\n",label_count);
      int* pix_counts =  new int[label_count];
      int* ind_counts =  new int[label_count];
      for(int i =0; i<label_count; i++){
         pix_counts[i] = IMG.SP(K+i]->N;
         ind_counts[i] = i;
            if(i>=num_SP) check_labels.addNode(K+i);
            %printf("SP: %d, %d pixels with label %d \n",K+i,IMG.SP(K+i]->N,label[IMG.SP(K+i]->pixels.getFirst()->getData()]);
      }
      U_quicksort<int>(pix_counts,ind_counts,0,label_count-1);

      %need to merge smallest one to its best neighbour
         %if(K==1719)
      %	mexErrMsgTxt("the redundant SP doesn't merge within the range...\n");

      double max_E;
      int max_k;
      arr(bool) neighbors = allocate_memory<bool>(K+label_count,false);
      for(int i =0; i<label_count-num_SP; i++)
      {
         max_k = -1;
         max_E = -DBL_MAX;
         move_merge_SP_propose_region( K+ind_counts[i], neighbors,check_labels,max_E,max_k);
         %printf("to merge %d,%d,%d\n",i,K+ind_counts[i],max_k);

         if(max_k==-1||max_E==-1)
            mexErrMsgTxt("the redundant SP doesn't merge within the range...\n");
         node = IMG.SP(K+ind_counts[i]]->pixels.getFirst();
         while(node !=NULL){
                  %if(K==689) printf("mmmerme:%d,%d\n",label[node->getData()],label[IMG.SP(max_k]->pixels.getFirst()->getData()]);
            IMG.label[node->getData()] = max_k;
            node = node->getNext();
         }
         IMG.SP(max_k]->merge_with(IMG.SP(K+ind_counts[i]],IMG.label, border_ptr, xdim, ydim);
         delete IMG.SP(K+ind_counts[i]];
         IMG.SP(K+ind_counts[i]] = NULL;
      }

      % relabel SPs
      %printf("relabel\n");
      max_k = K;
      for(int i =0; i<num_SP; i++)
      {
         while(IMG.SP(max_k]==NULL)
         {
            max_k++;
            if(max_k>N)
               mexErrMsgTxt("should be more nonempty SPs...\n");
         }
         if(max_k!=K+i)
         {
            %fetch stuff far behind to here
            IMG.SP(K+i] = IMG.SP(max_k];
            IMG.SP(max_k] = NULL;
            %printf("pair up %d,%d\n",(K+i),max_k);
            %relabel to the first two ...
            node = IMG.SP(K+i]->pixels.getFirst();
            while(node !=NULL)
            {
               IMG.label[node->getData()] = K+i;
               node = node->getNext();
            }
         }
         max_k++;
      }

      delete[] pix_counts;
      delete[] ind_counts;
      deallocate_memory(neighbors);
   }
   deallocate_memory(nlabels);
}



%