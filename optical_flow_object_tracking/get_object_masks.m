% returns a list of masked images
function object_masks = get_object_masks(params, root, flow_folder)
    addpath('util/');
    disp('Finding objects...');
    objects = find_objects([root flow_folder], params);
    if numel(objects)==0
        error('No objects found :(');
    end
    
    % combine similar objects
    disp('Combining objects...');
    objects = combine_objects(objects, params.sameness_threshold, params.time_inc);
    
    % interpolate objects over frames when they're lost and delete objects that don't persist
    disp('Interpolating objects...');
    objects = interpolate_objects(objects, params.min_frames_per_object, length(objects(1).found_masks));
    fprintf('%d objects found!\n', numel(objects));
    
    % fill in the masks between every [frame_inc]th mask
    object_masks = get_intermediate_masks(objects, params);
end       


% finds a mask for each image
function objects = find_objects(root, params)
    tolerance = 2*pi*params.percentage_tolerance/100;
    files = dir([root '*.mat']);
    max_f = numel(files);
            
    for f=1:max_f
        fprintf('Finding objects in flow file %d/%d\n', f, max_f);
        
        % get the boolean mask from the flow
        [boolean_mask, vx, vy] = get_boolean_flow_mask(f, root, params.cutoff_iters, files);
       
        %initialize empty struct array
        if f==1
            [xdim, ydim] = size(boolean_mask);
            new_obj.masks = false(xdim,ydim,max_f);
            new_obj.found_masks = false(max_f, 1);
            new_obj.theta_medians = zeros(max_f, 1);
            new_obj.dead = false;
            objects(1) = new_obj;
            objects(1) = [];
        end
        
        % make a separate mask for each object by finding connected areas in boolean_mask
        object_masks = find_connected_areas(boolean_mask);

        for l=1:size(object_masks, 3)
            
            % if this mask is below the minimum size, ignore it
            if sum(sum(object_masks(:,:,l))) < (xdim*ydim)*(params.min_object_percentage/100)
                continue;
            end
            
            % expand the mask to make sure we don't miss anything
            expanded_mask = expand_mask(object_masks(:,:,l), params.expansion);
            
            % find the median theta for this object
            [theta, theta_median, theta_offset] = get_theta_median(object_masks(:,:,l), vx, vy);

            % make things that point in the wrong direction false
            % make sure that we're accounting for the values wrapping around
            if theta_median-tolerance < -pi
                expanded_mask(theta>(theta_median+tolerance) & theta<(theta_median-tolerance+2*pi)) = false;
            elseif theta_median+tolerance > pi
                expanded_mask(theta>(theta_median+tolerance-2*pi) & theta<theta_median-tolerance) = false;
            else
                expanded_mask(theta>theta_median+tolerance | theta<theta_median-tolerance) = false;
            end
            
            % calculate the absolute median theta for this mask by removing the offset
            theta_median = get_absolute_theta(theta_median, theta_offset);
            
            % check to see if it matches another object
            max_sameness = 0; % initialize the max sameness
            max_sameness_object = 0;
            theta_used = theta_median; % default to the theta_median
            for o=1:numel(objects)
                % make sure the object isn't dead
                if objects(o).dead
                    continue;
                end
                
                % find the most recent found_mask in the given object
                last_idx = find(objects(o).found_masks(1:f-1), 1, 'last');
                if isempty(last_idx) % if this index doesn't exist, continue
                    continue;
                end
                
                % create a list of the last 5 thetas of o to check
                thetas_to_check = theta_median;
                ff = f-1;
                while ff > 0 && length(thetas_to_check) < 5
                    if objects(o).found_masks(ff)
                        thetas_to_check(end+1) = objects(o).theta_medians(ff);
                    end
                    ff = ff - 1;
                end

                % see which recent theta is best
                for t=1:length(thetas_to_check)
                    curr_sameness = sameness(expanded_mask, objects(o).masks(:,:,last_idx), thetas_to_check(t));
                    if curr_sameness > max_sameness
                        max_sameness = curr_sameness;
                        max_sameness_object = o;
                        theta_used = thetas_to_check(t); % store the theta that was used to set later
                    end
                end
            end

            % if the new object doesn't match any old objects
            if max_sameness < params.sameness_threshold
                objects(end+1) = new_obj; % create a new object
                max_sameness_object = numel(objects); % this is the object to add the current info to
                theta_used = theta_median;
            end

            % update the appropriate object with the new information
            objects(max_sameness_object).masks(:,:,f) = objects(max_sameness_object).masks(:,:,f) | expanded_mask;
            objects(max_sameness_object).theta_medians(f) = theta_used;
            objects(max_sameness_object).found_masks(f) = true;
        end
        
        % check to see if objects are "dead", i.e. haven't been updated for 2 seconds
        for o=1:numel(objects)
            if ~objects(o).dead && f - find(objects(o).found_masks, 1, 'last') >= 2 / params.time_inc
                objects(o).dead = true;
            end
        end
    end
end

% get a boolean mask showing the where the flow is highest, and also return
% the added flow from frames before and after
function [boolean_mask, final_vx, final_vy] = get_boolean_flow_mask(f, root, cutoff_iters, files)
    max_f = numel(files);
    for ff=f:2:min(max_f, f+2)
        if ff~=f
            old_boolean_mask = boolean_mask;
        end
        % flow from past
        load(fullfile(root, files(ff).name));
        if ff==f
            vx = -flow.bvx;
            vy = -flow.bvy;
        else
            final_vx = vx + flow.fvx;
            final_vy = vy + flow.fvy;
            vx = flow.fvx;
            vy = flow.fvy;
        end
        additive_flow = abs(vx) + abs(vy);

        %make a mask of the flow values
        for iter=1:cutoff_iters
            cutoff = mean(mean(additive_flow(additive_flow>0)));
            additive_flow(additive_flow<cutoff) = 0;
        end
        additive_flow(additive_flow>=cutoff) = 1;

        boolean_mask = false(size(additive_flow));
        boolean_mask(additive_flow==1) = true;
    end

    if f<=max_f-2
        boolean_mask = boolean_mask & old_boolean_mask;
    else
        final_vx = vx;
        final_vy = vy;
    end
end


% combine objects that are the same
function objects = combine_objects(objects, sameness_threshold, time_inc)
    o=1;
    while o <= numel(objects)
        max_sameness = 0;
        max_sameness_object = 0;
        theta_used = 0;
        theta_to_change = 0;
        for oo=o+1:numel(objects)
            last_oo = find(objects(oo).found_masks, 1);
            last_o = find(objects(o).found_masks(1:last_oo-1), 1, 'last');
            
            % make sure that the frames are not too far appart
            if isempty(last_o) || isempty(last_oo) || last_oo - last_o > 5 / time_inc
                continue;
            end
            
            % find the last few thetas from o and get the median of them
            thetas_to_check = objects(o).theta_medians(last_o);
            f = last_o-1;
            while f > 0 && length(thetas_to_check) < 4
                if objects(o).found_masks(f)
                    thetas_to_check(end+1) = objects(o).theta_medians(f);
                end
                f = f - 1;
            end
            [~, theta_median_o, theta_offset_o] = get_theta_median(thetas_to_check);
            theta_median_o = get_absolute_theta(theta_median_o, theta_offset_o);
            
            % find the last few medians from oo and get the median of them
            f = last_oo;
            oo_start = length(thetas_to_check)+1;
            while f <= length(objects(oo).found_masks) && length(thetas_to_check) < oo_start + 3
                if objects(oo).found_masks(f)
                    thetas_to_check(end+1) = objects(oo).theta_medians(f);
                end
                f = f + 1;
            end
            [~, theta_median_oo, theta_offset_oo] = get_theta_median(thetas_to_check(oo_start:end));
            theta_median_oo = get_absolute_theta(theta_median_oo, theta_offset_oo);

            % direction rating
            dir_rating = abs(theta_median_o - theta_median_oo);
            if dir_rating>pi
                dir_rating = 2*pi - dir_rating;
            end
            dir_rating = (pi - dir_rating) / pi;
            
            % check to see if we actually have a direction rating
            if isnan(dir_rating)
                continue;
            end

            % see which recent theta is best
            for t=1:length(thetas_to_check)
                % overlap rating
                curr_sameness = sameness(objects(oo).masks(:,:,last_oo), objects(o).masks(:,:,last_o), thetas_to_check(t)) * (dir_rating^2);
                if curr_sameness > max_sameness
                    max_sameness = curr_sameness;
                    max_sameness_object = oo;
                    theta_used = thetas_to_check(t);
                    theta_to_change = last_oo;
                end
            end
        end
        if max_sameness > sameness_threshold
            %merge the objects
            objects(o).found_masks = objects(o).found_masks | objects(max_sameness_object).found_masks;
            objects(o).theta_medians(objects(max_sameness_object).found_masks) = objects(max_sameness_object).theta_medians(objects(max_sameness_object).found_masks);
            objects(o).theta_medians(theta_to_change) = theta_used;
            objects(o).masks(:,:, objects(max_sameness_object).found_masks) = objects(max_sameness_object).masks(:, :, objects(max_sameness_object).found_masks);
            objects(max_sameness_object) = [];
        else
            o = o+1;
        end
    end
end


% interpolate objects into missing frames and delete objects that don't persist
function objects = interpolate_objects(objects, min_frames_per_object, max_f)
    o=1;
    while o <= numel(objects)
        % get rid of objects that don't persist
        if sum(objects(o).found_masks) < min_frames_per_object
            objects(o) = [];
        else
            % interpolate objects over frames where they're lost
            last = 0;
            f = find(objects(o).found_masks, 1);
            while f <= max_f
                while f <= max_f && objects(o).found_masks(f)
                    last = f;
                    f = f+1;
                end
                idx_ahead = find(objects(o).found_masks(f:end), 1);
                if idx_ahead
                    f_next = idx_ahead+f-1;
                    [~, row_offset, col_offset] = get_overlap(objects(o).theta_medians(f_next), objects(o).masks(:,:,last), objects(o).masks(:,:,f_next), false);
                    
                    for ff=1:idx_ahead-1
                        objects(o).masks(:,:,f) = get_offset_mask(objects(o).masks(:,:,last), round(row_offset * ff/idx_ahead), round(col_offset * ff/idx_ahead));
                        objects(o).theta_medians(f) = objects(o).theta_medians(f_next);
                        f = f+1;
                    end
                else
                    break;
                end
            end
            o=o+1;
        end
    end
end


% find the masks to fill in between every [frame_inc]th mask
function objects_final_masks = get_intermediate_masks(objects, params)
    %objects_final_masks = false(size(objects(1).masks,1), size(objects(1).masks,2), (size(objects(1).masks,3)+1)*params.frame_inc+1, numel(objects));
    objects_final_masks = false(size(objects(1).masks,1), size(objects(1).masks,2), size(objects(1).masks,3)*params.frame_inc, numel(objects));
    for o=1:numel(objects)
        final_masks = false(size(objects_final_masks, 1), size(objects_final_masks, 2), size(objects_final_masks, 3));
        for frame=1:size(objects(o).masks, 3)
            original_mask = objects(o).masks(:,:,frame);
            if frame~=1
                [before_overlap, row_offset_before, col_offset_before] = get_overlap(objects(o).theta_medians(frame), objects(o).masks(:,:,frame-1), original_mask, false);
                curr_mask = original_mask | before_overlap;
            end
            if frame < size(objects(o).masks,3)
                [after_overlap, row_offset_after, col_offset_after] = get_overlap(objects(o).theta_medians(frame+1), objects(o).masks(:,:,frame+1), original_mask, true);
                curr_mask = original_mask | after_overlap;
            end

            if frame~=1 && frame < size(objects(o).masks,3)
                curr_mask = (original_mask & after_overlap) | (original_mask & before_overlap) | (after_overlap & before_overlap);
            end
            
            % THIS IS THE CHANGE
            %curr_mask = original_mask;
            
            frames_back = params.frame_inc;
            frames_forward = 0;

            if frame==1
                row_offset = -row_offset_after;
                col_offset = -col_offset_after;
            elseif frame==size(objects(o).masks,3)
                row_offset = row_offset_before;
                col_offset = col_offset_before;
            else
                row_offset = mean([row_offset_before -row_offset_after]);
                col_offset = mean([col_offset_before -col_offset_after]);
            end

            for f=(-frames_back+1):frames_forward
                final_masks(:,:,frame*params.frame_inc+f) = get_offset_mask(curr_mask, round(f*row_offset/params.frame_inc), round(f*col_offset/params.frame_inc));
            end
        end
        objects_final_masks(:,:,:,o) = final_masks;
    end
end