% =============================================================================
% == split_move.cpp
% == --------------------------------------------------------------------------
% == A MEX interface to perform split moves on TSPs.
% == See m files for calling convention.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-05-2014
% =============================================================================

function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = move_split_SP(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_alive_dead_changed, IMG_N, IMG_new_SP, IMG_SP_changed, Nsp, index, forced)
    num_SP = 2;
    %fprintf('doin a split doe! on TSP %d\n', index);
    if ~(index > numel(IMG_SP) || isempty(IMG_SP(index).N) || IMG_SP(index).N == 0) && IMG_SP(index).N > num_SP
        new_ks = ones(num_SP, 1) * -1;
        ksplit = -1;
        max_E = -inf;

        %only work for num_SP == 2, so far
        [IMG_label, IMG_SP, max_E, ksplit, new_ks] = move_split_SP_propose(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_N, IMG_new_SP, index, num_SP, max_E, ksplit, new_ks);
        
        % update
        if max_E>0 || (forced && new_ks(1) > 0)
            %mexPrintf("split: %f,%d,%d\n",max_E,index,new_ks[0]);
            % update the labels first
            IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(new_ks(1)), index);

            % merge the super pixels
%            disp(IMG_SP(new_ks(1)).N);
            IMG_SP = U_merge_SPs(IMG_SP, IMG_label, index, new_ks(1)); %need to delete
            IMG_SP_changed(index) = true;
            
            if (ksplit<0)
                % splitting into a new super pixel
                IMG_SP_changed(new_ks(2)) = true;
                % move it to the right spot
                if (new_ks(2)~=IMG_K+1)
                    IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(new_ks(2)), IMG_K+1);
                    IMG_SP = U_merge_SPs(IMG_SP, IMG_label, IMG_K+1, new_ks(2));
                    IMG_SP_old(IMG_K+1) = IMG_SP_old(new_ks(2));
                end
                IMG_K = IMG_K + 1;
                IMG_max_UID = IMG_max_UID + 1;
            else
                % splitting into an old super pixel
                IMG_SP_changed(ksplit) = true;
                IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(new_ks(2)), ksplit);
                IMG_SP = U_merge_SPs(IMG_SP, IMG_label, ksplit, new_ks(2));
                IMG_alive_dead_changed = true;
            end
        elseif (ksplit~=-2)
            %Recover previous IMG_label
            IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(IMG_K+1), index);
            IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(IMG_K+2), index);
            IMG_SP_changed(index) = false;

            for i=1:num_SP
                IMG_SP = U_merge_SPs(IMG_SP, IMG_label, index, IMG_K+i);
            end
        else
            IMG_SP_changed(index) = false;
        end
    end
    IMG_SP = remove_excess_SPs(IMG_SP, IMG_K, Nsp);
end


function [IMG_label, IMG_SP, max_E, ksplit, new_ks] = move_split_SP_propose(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_N, IMG_new_SP, index, num_SP, max_E, ksplit, new_ks)
    [xdim, ydim] = size(IMG_label);
    num_iter = 5;
    SP_bbox = [xdim, 1, ydim, 1];

    % 1. Kmeans++
    [IMG_label, SP_bbox, broken] = U_Kmeans_plusplus(IMG_label, IMG_SP, IMG_data, IMG_SP_old, IMG_K, index, SP_bbox, num_SP, num_iter);
    
    if broken
        max_E = -inf;
        ksplit = -2;
        return;
    end

    % 2. create two new SP from the old one
    % label matrix is equivalent to SP.pixels
    % we will mess around label matrix for new proposal
    % and recover it from SP.pixels if it doesn't work out
    % merge small connected component in the label matrix
    [IMG_label, IMG_SP] = U_connect_newSP(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, IMG_N, IMG_SP_old, model_order_params, index, SP_bbox, num_SP);

    % new super pixels in K+1 and K+2... old super pixel in index
    IMG_SP(index) = SP_empty(IMG_SP(index));

%    if (option==-1)
    % (index, new);
    E = move_split_calc_delta(model_order_params, IMG_SP(index), IMG_new_SP, IMG_SP(IMG_K+1), IMG_SP(IMG_K+2));
    if (E>max_E)
        max_E = E;
        new_ks(1) = IMG_K+1;
        new_ks(2) = IMG_K+2;
        ksplit = -1;
    end
    % (new, index)
    E = move_split_calc_delta(model_order_params, IMG_SP(index), IMG_new_SP, IMG_SP(IMG_K+2), IMG_SP(IMG_K+1));
    if (E>max_E)
        max_E = E;
        new_ks(1) = IMG_K+2;
        new_ks(2) = IMG_K+1;
        ksplit = -1;
    end
    % (index, old_empty) && (old_empty, index)
    for ktest=1:IMG_K
        if (ktest~=index && ktest<=numel(IMG_SP) && IMG_SP(ktest).N==0)
            E = move_split_calc_delta(model_order_params, IMG_SP(index), IMG_SP(ktest), IMG_SP(IMG_K+1), IMG_SP(IMG_K+2));
            if (E>max_E)
                max_E = E;
                new_ks(1) = IMG_K+1;
                new_ks(2) = IMG_K+2;
                ksplit = ktest;
            end
            E = move_split_calc_delta(model_order_params, IMG_SP(index), IMG_SP(ktest), IMG_SP(IMG_K+2), IMG_SP(IMG_K+1));
            if (E>max_E)
                max_E = E;
                new_ks(1) = IMG_K+2;
                new_ks(2) = IMG_K+1;
                ksplit = ktest;
            end
        end
    end
end


function dist = U_dist(vec1, vec2)
   dist = sum((vec1 - vec2).^2);
   dist = dist/10000;
end


function [IMG_label, bbox, broken] = U_Kmeans_plusplus(IMG_label, IMG_SP, IMG_data, IMG_SP_old, IMG_K, index, bbox, num_SP, numiter)
    if (num_SP~=2)
        disp('Trying to split into more than 2');
    end
    
    broken = false;
    curr_SP = IMG_SP(index);
    true_pix = find(curr_SP.pixels)';
    xdim = size(IMG_label, 1);

    num_pix = curr_SP.N;
    distvec = zeros(length(curr_SP.pixels),1);
    klabels = zeros(num_pix,1);
    center = zeros(num_SP, 5);

    %1. kmeans ++ initialization
    %first cener pt

    for tmp_pos_index=1:length(true_pix)
        tmp_pos = true_pix(tmp_pos_index);
        % get the Bounding Box of curr_SP
        [x, y] = get_x_and_y_from_index(tmp_pos, xdim);
        bbox(1) = min(x, bbox(1));
        bbox(2) = max(x, bbox(2));
        bbox(3) = min(y, bbox(3));
        bbox(4) = max(y, bbox(4));
    end

    old_split1 = IMG_SP_old(index);
    if ~old_split1 % new
        % first is new
        seed = true_pix(randi(length(true_pix)));

        [x, y] = get_x_and_y_from_index(seed, xdim);
        if (IMG_label(x, y)~=index)
            disp('inconsistency about cluster label');
        end

        center(1,:) = IMG_data(seed,:);
        
        for pix_index=1:length(true_pix)
            pix = true_pix(pix_index);
            distvec(pix) = U_dist(IMG_data(pix,:), center(1,:));
        end

        [~, seed2] = max(distvec);
        [x, y] = get_x_and_y_from_index(seed2, xdim);
        if (IMG_label(x, y)~=index)
            disp('inconsistency about cluster label');
        end
        center(2,:) = IMG_data(seed2,:);

    else %old
        center(1,1:2) = NormalD_calc_mean(curr_SP.pos);
        center(1,3:5) = NormalD_calc_mean(curr_SP.app);

        for pix_index=1:length(true_pix)
            pix = true_pix(pix_index);
            distvec(pix) = U_dist(IMG_data(pix,:), center(1,:));
        end

        % second is new
        [~, seed2] = max(distvec);
        [x, y] = get_x_and_y_from_index(seed2, xdim);
        if (IMG_label(x, y)~=index)
            disp('inconsistency about cluster label');
        end
        center(2,:) = IMG_data(seed2,:);
     end

    %2. kmeans ++ iterations
    change = false;
    for itr=1:numiter
        distvec = inf(length(curr_SP.pixels),1);
        for pix_index=1:length(true_pix)
            pix = true_pix(pix_index);
            for i=1:num_SP
                tmp_dist = U_dist(IMG_data(pix,:), center(i,:));
                if(tmp_dist < distvec(pix) )
                    distvec(pix) = tmp_dist;
                    klabels(pix) = i;
                    change = true;
                end
            end
        end

        if ~change
            %no change happened... Kmeans totally stuck
            break;
        end

        center = zeros(num_SP, 5);
        SP_sz = zeros(num_SP, 1);


        for pix_index=1:length(true_pix)
            pix = true_pix(pix_index);

            % klabels[n]==0 && old_split1 then don't update
            % klabels[n]==1 && old_split2 then don't update
            % if (klabels[n]!=old_split1-1 && klabels[n]!=old_split2-1)
            if (klabels(pix)~=1 || ~old_split1)
                center(klabels(pix), :) = center(klabels(pix), :) + IMG_data(pix, :);
            end
            SP_sz(klabels(pix)) = SP_sz(klabels(pix))+1;
        end

        for k=1:num_SP
            if SP_sz(k)>0
                if (k==1 && ~old_split1) || k==2
                    center(k, :) = center(k, :) / SP_sz(k);
                end
            else
                if old_split1
                    broken = true;
                else
                    disp('one cluster removed... shouldnt happen');
                end
            end
        end
    end

    if ~broken
        % change label accordingly
        for pix_index=1:length(true_pix)
            pix = true_pix(pix_index);
            [x, y] = get_x_and_y_from_index(pix, xdim);
            IMG_label(x, y) = IMG_K+klabels(pix);
        end
    end
end


% KMerge Version 2: grow region

function [IMG_label, IMG_SP] = U_connect_newSP(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, IMG_N, IMG_SP_old, model_order_params, old_k, bbox, num_SP)
    min_x = bbox(1);
    min_y = bbox(3);
    max_x = bbox(2);
    max_y = bbox(4);
    dx4 = [-1,  0,  1,  0];
    dy4 = [ 0, -1,  0,  1];

    [xdim, ydim] = size(IMG_label);
    new_label = false(xdim, ydim);
    label_count = 0;

    check_labels = (1:num_SP) + IMG_K;

    for x=min_x:max_x
        for y=min_y:max_y
            tmp_ind = get_index_from_x_and_y(x, y, xdim);
            curLabel = IMG_label(x, y);
            if any(curLabel==check_labels) && ~new_label(x, y)
                label_count = label_count + 1;
                k = IMG_K + label_count;
                new_label(x, y) = true;
                
                if ~(k > numel(IMG_SP) || isempty(IMG_SP(k).N) || IMG_SP(k).N == 0)
                    disp('SP should be null..');
                end
                IMG_SP(k) = new_SP(IMG_new_pos, IMG_new_app, IMG_max_UID, [0, 0], IMG_N, IMG_max_SPs);
                
                %can't decide whether it will be border yet
                pixel_bfs=zeros(size(IMG_SP(k).pixels));
                pixel_bfs_write=1;
                pixel_bfs_read=1;
                
                IMG_SP(k) = SP_add_pixel(IMG_SP(k), IMG_data, tmp_ind, false, IMG_boundary_mask(x, y));
                pixel_bfs(pixel_bfs_write) = tmp_ind;
                pixel_bfs_write = pixel_bfs_write+1;
                
                while(pixel_bfs_read < pixel_bfs_write)
                    pix = pixel_bfs(pixel_bfs_read);
                    [pix_x, pix_y] = get_x_and_y_from_index(pix, xdim);
                    for n=1:4
                        new_x = pix_x+dx4(n);
                        new_y = pix_y+dy4(n);
                        if (new_x >= min_x && new_x <= max_x) && (new_y >= min_y && new_y <= max_y)
                            new_ind = get_index_from_x_and_y(new_x, new_y, xdim);
                            if ~new_label(new_x, new_y) && IMG_label(new_x, new_y) == curLabel
                                % should always update labels before adding pixel, otherwise
                                % U_check_border_pix will be wrong!
                                new_label(new_x, new_y) = true;
                                IMG_label(new_x, new_y) = k;
                                IMG_SP(k) = SP_add_pixel(IMG_SP(k), IMG_data, new_ind, false, IMG_boundary_mask(new_x, new_y));
                                pixel_bfs(pixel_bfs_write) = new_ind;
                                pixel_bfs_write = pixel_bfs_write+1;
                            end
                        end
                    end
                    pixel_bfs_read = pixel_bfs_read+1;
                end
                IMG_label(x, y) = k;
            end
        end
    end

    % Now is the time to clean up the border pixels and neighbor ids
    
    % remove all neighbors of the original index
    for k=1:IMG_K
        IMG_SP(k).neighbors(old_k) = 0;
    end
    
    
    % loop through every pixel in the bounding box plus one pixel all around
    % set border and neighbor values
    for x=max(1, min_x-1):min(xdim, max_x+1)
        for y=max(1, min_y-1):min(ydim, max_y+1)
            k = IMG_label(x, y);
            if k>0
                index = get_index_from_x_and_y(x, y, xdim);
                is_border = U_check_border_pix(IMG_label, index, k);
                IMG_SP(k).borders(index) = is_border;
                if is_border
                    % 1:IMG_K are already set up for neighbor 1:IMG_K, so
                    % only update them for the new SPs
                    if k<=IMG_K
                        old_neighbors = IMG_SP(k).neighbors(1:IMG_K);
                    end
                    IMG_SP(k) = SP_update_neighbors_add_self(IMG_SP(k), IMG_label, index);
                    if k<=IMG_K
                        IMG_SP(k).neighbors(1:IMG_K) = old_neighbors;
                    end
                end
            end
        end
    end
    
    if (label_count>num_SP)
        pix_counts = zeros(label_count, 1);
        for i=1:label_count
            pix_counts(i) = IMG_SP(IMG_K+i).N;
        end
        check_labels = (IMG_K+1):(IMG_K+label_count);
        
        [~, ind_counts] = sort(pix_counts, 1, 'ascend');
        
        neighbors = false(IMG_K+label_count, 1);
        for i=1:label_count-num_SP
            max_k = -1;
            max_E = -inf;
            
            [max_E, max_k] = move_merge_SP_propose_region(IMG_label, IMG_SP, IMG_SP_old, model_order_params, IMG_K+ind_counts(i), neighbors, check_labels, max_E, max_k);
            
            if(max_k==-1 || max_E==-1)
                disp('the redundant SP doesnt merge within the range');
                save('connect_newSP.mat');
            end
            
            IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(IMG_K+ind_counts(i)), max_k);
            IMG_SP = U_merge_SPs(IMG_SP, IMG_label, max_k, IMG_K+ind_counts(i));
        end

        % relabel SPs
        %printf("relabel\n");
        max_k = IMG_K+1;
        for i=1:num_SP
            while (max_k > numel(IMG_SP) || isempty(IMG_SP(max_k).N) || IMG_SP(max_k).N == 0)
                max_k = max_k+1;
                if(max_k>IMG_max_SPs)
                    disp('there should be more nonempty SPs..\n');
                end
            end
            if max_k~=IMG_K+i
                %fetch stuff far behind to here
                IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(max_k), IMG_K+i);
                IMG_SP = U_merge_SPs(IMG_SP, IMG_label, IMG_K+i, max_k);
                %printf("pair up %d,%d\n",(K+i),max_k);
                %relabel to the first two ...
            end
            max_k = max_k+1;
        end
    end
end


function [max_E, max_k] = move_merge_SP_propose_region(IMG_label, IMG_SP, IMG_SP_old, model_order_params, k, neighbors, check_labels, max_E, max_k)

    neighbors = U_find_border_SP(IMG_label, IMG_SP, k, neighbors);

    % loop through all neighbors
    neighbor_intersection = intersect(find(neighbors)', check_labels);
    for merge_k_index=1:length(neighbor_intersection);
        merge_k = neighbor_intersection(merge_k_index);
        % calculate the energy
        tmp_E = U_move_merge_calc_delta(IMG_SP(k), IMG_SP(merge_k), IMG_SP_old(k), IMG_SP_old(merge_k), model_order_params);
        if(tmp_E>max_E)
            max_E = tmp_E;
            max_k = merge_k;
        end
    end
            
    if (max_k==-1)
        disp('move_merge_SP_propose_region: No neighbour found');
        save('propose.mat');
    end
end

function IMG_SP = remove_excess_SPs(IMG_SP, IMG_K, Nsp)
    %REMOVE EXCESS SPs
    curr_SP = numel(IMG_SP);
    while (isempty(IMG_SP(curr_SP).N) || IMG_SP(curr_SP).N == 0) && curr_SP>Nsp
        IMG_SP(curr_SP) = [];
        curr_SP = curr_SP - 1;
    end
    if IMG_K~=numel(IMG_SP)
        disp('it done broke');
    end
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
function prob = move_split_calc_delta(model_order_params, SP1, SP2, new_SP1, new_SP2)
    prob = SP_log_likelihood_test_merge1(SP1, new_SP1);
    prob = prob + SP_log_likelihood_test_merge1(SP2, new_SP2);
    prob = prob - SP_log_likelihood_test_merge2(SP1, new_SP1, new_SP2);
    prob = prob - SP2.log_likelihood;

    % split
    size1 = SP1.N + new_SP1.N;
    size_var1 = (model_order_params.area - size1/2)*size1/model_order_params.area_var;
    size2 = SP2.N + new_SP2.N;
    size_var2 = (model_order_params.area - size2/2)*size2/model_order_params.area_var;

    % not split
    size3 = SP1.N + new_SP1.N + new_SP2.N;
    size_var3 = (model_order_params.area - size3/2)*size3/model_order_params.area_var;
    size4 = SP2.N;
    size_var4 = (model_order_params.area - size4/2)*size4/model_order_params.area_var;

    prob = prob + size_var1 + size_var2 - size_var3 - size_var4;
end