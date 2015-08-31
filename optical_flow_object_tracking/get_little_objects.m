clear;
addpath('util/');

load('object_masks.mat');

root = 'sequences/hop/';
flow_folder = 'little_flow/';
img_suffix = 'jpg';

% compute the optical flows for the frames
compute_flows(root, img_suffix, flow_folder, 1, true);

% get flow files
files = dir([strcat(root, flow_folder) '*.mat']);

params.expansion = 40;
params.min_frames_per_object = 10;
params.sameness_threshold = 0.2;

%initialize objects
[xdim, ydim] = size(object_masks(:,:,1,1));
max_f = size(object_masks, 3);
new_obj.masks = false(xdim,ydim,max_f);
new_obj.found_masks = false(max_f, 1);
new_obj.dead = false;

for o=1:size(object_masks, 4)
    centers = zeros(size(object_masks, 3), 2);
    for c = 1:size(centers, 1)
        centers(c,:) = get_center(object_masks(:,:,c,o));
    end
    clearvars little_objects;
    little_objects(1) = new_obj;
    little_objects(1) = [];
    for f = 1:max_f
        frame_index = f+1+params.frame_inc;
        fprintf('object %d, frame %d\n', o, frame_index);
        
        orig_mask = object_masks(:,:,f,o);
        if sum(orig_mask(:)) == 0
            continue;
        end
        
        for ff=frame_index-1:frame_index
            % flow from past
            %load(strcat(root, flow_folder, 'frame', num2str(ff, '%04d'), '_flow.mat'));
            load(fullfile(strcat(root, flow_folder), files(ff).name));
            if ff~=frame_index
                vx = -flow.bvx;
                vy = -flow.bvy;
            else
                old_additive_flow = additive_flow;
                final_vx = vx + flow.fvx;
                final_vy = vy + flow.fvy;
                vx = flow.fvx;
                vy = flow.fvy;
            end
            additive_flow = abs(vx) + abs(vy);
        end

        additive_flow = additive_flow + old_additive_flow;
        additive_flow_save = additive_flow;
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
%             max_sameness = 0; % initialize the max sameness
%             max_sameness_object = 0;
            found = false;
            oo = 1;
            first_object = 0;
            while oo <= numel(little_objects)
%                fprintf('oo = %d, numel(little_objects) = %d\n', oo, numel(little_objects));
                % make sure the object isn't dead
                if little_objects(oo).dead
                    oo = oo+1;
                    continue;
                end
                
                % find the most recent found_mask in the given object
                last_idx = find(little_objects(oo).found_masks(1:f-1), 1, 'last');
                if isempty(last_idx) % if this index doesn't exist, continue
                    oo = oo+1;
                    continue;
                end

                offset = centers(f,:) - centers(last_idx,:);
                calculated_theta = atan2(offset(1), offset(2));
                
                % see which recent theta is best
                %curr_sameness = sameness(curr_mask, little_objects(oo).masks(:,:,last_idx), calculated_theta);
                offset_mask = get_offset_mask(little_objects(oo).masks(:,:,last_idx), offset(1), offset(2));
                curr_sameness = sum( sum( curr_mask & offset_mask ) ) / min(sum( curr_mask(:) ), sum( sum(little_objects(oo).masks(:,:,last_idx)) ));
%                 if curr_sameness > max_sameness
%                     max_sameness = curr_sameness;
%                     max_sameness_object = oo;
%                 end
                if curr_sameness > params.sameness_threshold
                    found = true;
                    if first_object == 0
                        first_object = oo;
                        little_objects(oo).masks(:,:,f) = little_objects(oo).masks(:,:,f) | curr_mask;
                        little_objects(oo).found_masks(f) = true;
                    else
                        little_objects(first_object).found_masks = little_objects(first_object).found_masks | little_objects(oo).found_masks;
                        little_objects(first_object).masks = little_objects(first_object).masks | little_objects(oo).masks;
                        little_objects(oo) = []; % delete unecessary object
                        oo = oo-1;
                    end
                end
                oo = oo+1;
            end

            if ~found
                little_objects(end+1) = new_obj;
                little_objects(end).masks(:,:,f) = curr_mask;
                little_objects(end).found_masks(f) = true;
            end
            % if the new object doesn't match any old objects
%             if max_sameness < params.sameness_threshold
%                 little_objects(end+1) = new_obj; % create a new object
%                 max_sameness_object = numel(little_objects); % this is the object to add the current info to
%             end

            % update the appropriate object with the new information
%             little_objects(max_sameness_object).masks(:,:,f) = curr_mask;
%             little_objects(max_sameness_object).found_masks(f) = true;
        end
        
        % check to see if objects are "dead", i.e. haven't been updated for 2 seconds
        for oo=1:numel(little_objects)
            if ~little_objects(oo).dead && f - find(little_objects(oo).found_masks, 1, 'last') >= 3 / 0.1
                little_objects(oo).dead = true;
            end
        end
    end
    
    %combine objects
%     oo = 1;
%     while oo <= numel(little_objects) % loop through each object
%         ooo = 1;
%         while ooo <= numel(little_objects) % loop through each other object
%             if ooo == oo
%                 continue; % don't check against itself
%             end
%             % need to check if any have the same mask in the same frame,
%             % not just a mask in the same frame
%             if any(little_objects(oo).found_masks & little_objects(ooo).found_masks)
%                 % if yes, combine found_masks and masks
%                 little_objects(oo).found_masks = little_objects(ooo).found_masks | little_objects(oo).found_masks;
%                 little_objects(oo).masks = little_objects(ooo).masks | little_objects(oo).masks;
%                 
%                 little_objects(ooo) = []; % delete unecessary object
%                 
%                 if ooo < oo % take care of indices
%                     oo = oo - 1;
%                 end
%             else
%                 ooo = ooo+1;
%             end
%         end
%         oo = oo+1;
%     end
        
    
    % get rid of non-persistent objects
    % add the little objects to the big object
    oo = 1;
    while oo <= numel(little_objects)
        if sum(little_objects(oo).found_masks) < params.min_frames_per_object
            little_objects(oo) = [];
        else
            object_masks(:,:,:,o) = object_masks(:,:,:,o) | little_objects(oo).masks;
            oo = oo + 1;
        end
    end
end
create_videos(object_masks, 'sequences/wheelock3/', params.frame_inc, 'results/wheelock3_new.avi', 10, false);

little_object_masks = false(size(object_masks, 1), size(object_masks, 2), size(object_masks, 3), numel(little_objects));
for o=1:numel(little_objects)
    little_object_masks(:,:,:,o) = little_objects(o).masks;
end
save('little_objects.mat', 'object_masks', 'little_object_masks', 'little_objects', 'params');
create_videos(little_object_masks, 'sequences/wheelock3/', params.frame_inc, 'results/wheelock3_little.avi', 10, true);