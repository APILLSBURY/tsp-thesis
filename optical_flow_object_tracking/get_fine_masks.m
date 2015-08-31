function object_masks = get_fine_masks(object_masks, params, root, flow_folder)

    % get flow files
    files = dir([strcat(root, flow_folder) '*.mat']);

    %initialize objects
    [xdim, ydim] = size(object_masks(:,:,1,1));
    max_f = size(object_masks, 3);
    new_obj.masks = false(xdim,ydim,max_f);
    new_obj.found_masks = false(max_f, 1);
    new_obj.dead = false;

    for o=1:size(object_masks, 4)
        centers = zeros(size(object_masks, 3), 2);
        
        % find the centers of all the masks
        for c = 1:size(centers, 1)
            centers(c,:) = get_center(object_masks(:,:,c,o));
        end
        
        clearvars fine_objects;
        fine_objects(1) = new_obj;
        fine_objects(1) = [];
        for f = 1:max_f
            frame_index = f+1+params.frame_inc;
            fprintf('object %d, frame %d\n', o, frame_index);

            orig_mask = object_masks(:,:,f,o);
            if sum(orig_mask(:)) == 0
                continue;
            end

            for ff=frame_index-1:frame_index
                % flow from past
                load(fullfile(strcat(root, flow_folder), files(ff).name));
                if ff~=frame_index
                    vx = -flow.bvx;
                    vy = -flow.bvy;
                else
                    old_additive_flow = additive_flow;
                    vx = flow.fvx;
                    vy = flow.fvy;
                end
                additive_flow = abs(vx) + abs(vy);
            end

            additive_flow = additive_flow + old_additive_flow;
            surrounding_mask = expand_mask(orig_mask, params.expansion) & ~orig_mask;
            additive_flow(~surrounding_mask) = 0;

            surrounding_area = sum(surrounding_mask(:));

            %make a mask of the flow values
            for iter=1:params.cutoff_iters
                cutoff = mean(mean(additive_flow(additive_flow>0)));
                additive_flow(additive_flow<cutoff) = 0;
            end

            boolean_mask = false(size(additive_flow));
            boolean_mask(additive_flow>0) = true;

            possible_object_masks = find_connected_areas(boolean_mask);

            for l=1:size(possible_object_masks, 3)
                curr_mask = possible_object_masks(:,:,l);
                if sum(curr_mask(:)) < params.min_object_percentage/100 * surrounding_area
                    continue;
                end

                % find the median theta for this object
                [theta, theta_median, ~] = get_theta_median(curr_mask, vx, vy);

                % make things that point in the wrong direction false
                % make sure that we're accounting for the values wrapping around
                tolerance = 2*pi*params.percentage_tolerance/100;
                if theta_median-tolerance < -pi
                    curr_mask(theta>(theta_median+tolerance) & theta<(theta_median-tolerance+2*pi)) = false;
                elseif theta_median+tolerance > pi
                    curr_mask(theta>(theta_median+tolerance-2*pi) & theta<theta_median-tolerance) = false;
                else
                    curr_mask(theta>theta_median+tolerance | theta<theta_median-tolerance) = false;
                end

                % check to see if it matches another object
                found = false;
                oo = 1;
                first_object = 0;
                while oo <= numel(fine_objects)
                    % make sure the object isn't dead
                    if fine_objects(oo).dead
                        oo = oo+1;
                        continue;
                    end

                    % find the most recent found_mask in the given object
                    last_idx = find(fine_objects(oo).found_masks(1:f-1), 1, 'last');
                    if isempty(last_idx) % if this index doesn't exist, continue
                        oo = oo+1;
                        continue;
                    end


                    % get the offset mask by taking translating it the
                    % difference between the object centers, then use that
                    % to find the sameness
                    offset = centers(f,:) - centers(last_idx,:);
                    offset_mask = get_offset_mask(fine_objects(oo).masks(:,:,last_idx), offset(1), offset(2));
                    curr_sameness = sum( sum( curr_mask & offset_mask ) ) / min(sum( curr_mask(:) ), sum( sum(fine_objects(oo).masks(:,:,last_idx)) ));

                    % if our mask is similar enough to be part of this
                    % object, add it to this object
                    if curr_sameness > params.sameness_threshold
                        found = true;
                        % if this is the first object we're adding this mask to, just add it
                        if first_object == 0
                            first_object = oo;
                            fine_objects(oo).masks(:,:,f) = fine_objects(oo).masks(:,:,f) | curr_mask;
                            fine_objects(oo).found_masks(f) = true;
                        % if this is not the first object we're adding this
                        % mask to, combine this object with the first one
                        % and add the mask to the combined object
                        else
                            fine_objects(first_object).found_masks = fine_objects(first_object).found_masks | fine_objects(oo).found_masks;
                            fine_objects(first_object).masks = fine_objects(first_object).masks | fine_objects(oo).masks;
                            fine_objects(oo) = []; % delete unecessary object
                            oo = oo-1;
                        end
                    end
                    oo = oo+1;
                end

                % if it's not the same as any previously found objects,
                % create a new object for it
                if ~found
                    fine_objects(end+1) = new_obj;
                    fine_objects(end).masks(:,:,f) = curr_mask;
                    fine_objects(end).found_masks(f) = true;
                end
            end

            % check to see if objects are "dead", i.e. haven't been updated for 2 seconds
            for oo=1:numel(fine_objects)
                if ~fine_objects(oo).dead && f - find(fine_objects(oo).found_masks, 1, 'last') >= 3 / 0.1
                    fine_objects(oo).dead = true;
                end
            end
        end

        % get rid of non-persistent objects
        % add the fine objects to the big object
        oo = 1;
        while oo <= numel(fine_objects)
            if sum(fine_objects(oo).found_masks) < params.min_frames_per_object
                fine_objects(oo) = [];
            else
                object_masks(:,:,:,o) = object_masks(:,:,:,o) | fine_objects(oo).masks;
                oo = oo + 1;
            end
        end
    end
%     % save just the fine masks so we can make a video of them
%     fine_object_masks = false(size(object_masks, 1), size(object_masks, 2), size(object_masks, 3), numel(fine_objects));
%     for o=1:numel(fine_objects)
%         fine_object_masks(:,:,:,o) = fine_objects(o).masks;
%     end
%     save('fine_objects.mat', 'object_masks', 'fine_object_masks', 'fine_objects', 'params');
%     create_videos(little_object_masks, root, params.frame_inc, 'results/wheelock3_little.avi', 10, true);
end