% returns a list of masked images
function [ims] = get_masked_images(frame_inc, root, flow_folder, expansion, percentage_tolerance, min_object_percentage, cutoff_iters, test)
    disp('Getting masked images...');
    [masks, theta_medians] = find_masks([root flow_folder], expansion, percentage_tolerance, min_object_percentage, cutoff_iters);
    if test
        final_masks = masks;
    else
        if size(masks, 3) < 2
            error('not enough masks');
        end
        final_masks = false(size(masks,1), size(masks,2), (size(masks,3)+1)*frame_inc);
        for frame=1:size(masks, 3)
            
            curr_mask = masks(:,:,frame);
            if frame~=1
                [before_overlap, row_offset_before, col_offset_before] = get_overlap(theta_medians(frame), masks(:,:,frame-1), masks(:,:,frame), true);
                curr_mask = curr_mask | before_overlap;
            end
            if frame < size(masks,3)
                [after_overlap, row_offset_after, col_offset_after] = get_overlap(theta_medians(frame+1), masks(:,:,frame+1), masks(:,:,frame), false);
                curr_mask = curr_mask | after_overlap;
            end

            frames_back = floor((frame_inc-1)/2);
            frames_forward = ceil((frame_inc-1)/2);
            if frame==1
                frames_back = frame_inc;
                row_offset = -row_offset_after;
                col_offset = -col_offset_after;
            elseif frame==size(masks,3)
                frames_forward = frame_inc;
                row_offset = row_offset_before;
                col_offset = col_offset_before;
            else
                row_offset = mean([row_offset_before -row_offset_after]);
                col_offset = mean([col_offset_before -col_offset_after]);
            end

            for f=-frames_back:frames_forward
                final_masks(:,:,frame*frame_inc+f+1) = get_offset_mask(curr_mask, round(f*row_offset/frame_inc), round(f*col_offset/frame_inc));
            end
        end
    end
    ims = mask_images(final_masks, root, frame_inc);
end


% takes a list of masks and applies them to the images
function ims = mask_images(final_masks, root, frame_inc)
    files = dir([root '*.jpg']);
    ims = uint8(zeros(size(final_masks, 1), size(final_masks, 2), 3, min(numel(files)-frame_inc, size(final_masks, 3))));
    for f=(1+frame_inc):(size(ims,4)+frame_inc)
        findex = f-frame_inc;
        ims(:,:,:,findex) = imread(fullfile(root,files(f).name));
        mask = double(final_masks(:,:,findex));
        mask(mask==0) = 0.15;
        
        %apply the mask to the image
        ims(:,:,1,findex) = uint8(double(ims(:,:,1,findex)) .* mask);
        ims(:,:,2,findex) = uint8(double(ims(:,:,2,findex)) .* mask);
        ims(:,:,3,findex) = uint8(double(ims(:,:,3,findex)) .* mask);
    end
end


% finds the maximum overlap between a moving mask and a stationary mask
% given the direction that the moving mask is going
% USAGE: theta_median should be the theta for whichever is the later mask
% set reverse_theta to true when stationary mask is the later mask
function [max_offset_mask, max_row_offset, max_col_offset] = get_overlap(theta_median, moving_mask, stationary_mask, reverse_theta)
    if reverse_theta
        theta_median = theta_median + pi;
        if theta_median>pi
            theta_median = theta_median - 2*pi;
        end
    end

    [xdim, ydim] = size(moving_mask);
    max_overlap = -1;
    max_offset_mask = false(xdim, ydim);
    max_row_offset = 0;
    max_col_offset = 0;
    rho = 0;
    found_overlap = false;
    while rho==0 || sum(offset_mask(:))>0 && ~(overlap==0 && found_overlap)
        [col_offset, row_offset] = pol2cart(theta_median, rho);
        row_offset = round(row_offset);
        col_offset = round(col_offset);
        offset_mask = get_offset_mask(moving_mask, row_offset, col_offset);
        combined_mask = offset_mask & stationary_mask;
        overlap = sum(combined_mask(:));
        if overlap>0
            found_overlap = true;
        end
        if overlap>max_overlap
            max_overlap = overlap;
            max_offset_mask = offset_mask;
            max_row_offset = row_offset;
            max_col_offset = col_offset;
        end
        rho = rho + 1;
    end
end


% finds the mask that is offset by a given number of rows and columns
function offset_mask = get_offset_mask(moving_mask, row_offset, col_offset)
    offset_mask = false(size(moving_mask));
    offset_mask(max(1+row_offset, 1):min(end, end+row_offset), max(1+col_offset, 1):min(end, end+col_offset)) = moving_mask(max(1, 1-row_offset):min(end-row_offset, end), max(1, 1-col_offset):min(end-col_offset, end));
end       


% finds a mask for each image
function [masks, theta_medians] = find_masks(root, expansion, percentage_tolerance, min_object_percentage, cutoff_iters)

    curr_sameness = 0; % need to have it defined for f==1
    
    tolerance = 2*pi*percentage_tolerance/100;
   
    files = dir([root '*.mat']);
    masks = [];
    theta_medians = zeros(numel(files));
    
    max_f = numel(files);
    
    for f=1:max_f
        max_sameness = 0; % initialize the max sameness
        
        % flow from past
        load(fullfile(root, files(f).name));
        vx1 = flow.bvx;
        vy1 = flow.bvy;
        [xdim, ydim] = size(vx1);
        additive_flow = abs(vx1)+abs(vy1);

        % create the masks matrix
        if f==1
            masks = zeros(xdim, ydim, max_f);
        end
        
        %make a mask of the flow values
        additive_flow = additive_flow ./ max(additive_flow(:));
        for iter=1:cutoff_iters
            cutoff = mean(mean(additive_flow(additive_flow>0)));
            additive_flow(additive_flow<cutoff) = 0;
        end
        additive_flow(additive_flow>=cutoff) = 1;

        boolean_mask = false(size(additive_flow));
        boolean_mask(additive_flow==1) = true;

        % make a separate mask for each object by finding connected areas
        % in boolean_mask
        object_masks = find_levels(boolean_mask);

        for l=1:size(object_masks, 3)
            % if this mask is below the minimum size, ignore it
            if sum(sum(object_masks(:,:,l))) < (xdim*ydim)*(min_object_percentage/100)
                continue;
            end
            
            % expand the mask to make sure we don't miss anything
            expanded_mask = expand_mask(object_masks(:,:,l), expansion);
            
            % find the median theta for this object
            [theta, theta_median, theta_offset] = get_theta_median(object_masks(:,:,l), vx1, vy1);

            % make things that point in the wrong direction false
            % make sure that we're accounting for the values wrapping around
            if theta_median < tolerance-pi
                expanded_mask(theta>theta_median+tolerance & theta<theta_median-tolerance+2*pi) = false;
            elseif theta_median > pi-tolerance
                expanded_mask(theta>theta_median+tolerance-2*pi & theta<theta_median-tolerance) = false;
            else
                expanded_mask(theta>theta_median+tolerance | theta<theta_median-tolerance) = false;
            end
            
            % calculate the absolute median theta for this mask by removing the offset
            theta_median = theta_median - theta_offset;
            if theta_median>pi
                theta_median = theta_median - 2*pi;
            end
            if theta_median<-pi
                theta_median = theta_median + 2*pi;
            end
            
            % check to see if it's actually an object (i.e. not just
            % scattered dots)
            % find bounding box
            %   find index of first and last row and column
            %[row, col, ~] = find(expanded_mask);
            %if sum(expanded_mask(:)) > (1 + max(row) - min(row)) * (1 + max(col) - min(col)) / 4
            
            col_std = std(expanded_mask, 0, 1);
            row_std = std(expanded_mask, 0, 2);

            if mean([col_std(col_std>0) row_std(row_std>0)']) > 0.23
                imagesc(expanded_mask);
                % check to see if this is the right mask
                if f==1 || sum(sum(masks(:,:,f-1)))==0
                    right_mask = sum(expanded_mask(:)) > sum(sum(masks(:,:,f)));
                else
                    curr_sameness = sameness(expanded_mask, masks(:,:,f-1), theta_median, theta_medians(f-1));
                    right_mask = curr_sameness > max_sameness;
                end
                if right_mask
                    max_sameness = curr_sameness;
                    theta_medians(f) = theta_median;
                    masks(:,:,f) = expanded_mask;
                end
            else
                imagesc(expanded_mask);
                fprintf('');
            end
        end
        
        %if we're below a sameness threshold, just use the previous mask
        if f>1 && max_sameness<0.2 && sum(sum(masks(:,:,f-1))) > 0
            % if there's nothing like the first frame, assume that there was nothing good in the first frame
            if f==2
                theta_medians(1) = 0;
                masks(:,:,1) = false(xdim, ydim);
            end
            theta_medians(f) = theta_medians(f-1);
            masks(:,:,f) = masks(:,:,f-1);
        end
    end
end

% takes a boolean mask and implements a depth first search to find all
% connected areas
function object_mask = find_levels(boolean_mask)
    [xdim, ydim] = size(boolean_mask);
    object_mask = logical.empty(xdim, ydim, 0);
    x=1;
    level = 0;
    while x<=xdim
        temp_y = find(boolean_mask(x,:), 1);
        if isempty(temp_y)
            x = x+1;
            continue;
        end
        temp_x = x;
        bfs_write = 1;
        bfs_read = 1;
        bfs = zeros(xdim*ydim, 2);
        bfs(bfs_write, :) = [temp_x, temp_y];
        level = level + 1;
        object_mask(:, :, level) = false(xdim, ydim); % create a new level
        object_mask(temp_x, temp_y, level) = true;
        boolean_mask(temp_x, temp_y) = false;
        bfs_write = bfs_write+1;
        while bfs_read < bfs_write
            temp_x = bfs(bfs_read, 1);
            temp_y = bfs(bfs_read, 2);
            if (temp_x>1 && boolean_mask(temp_x-1, temp_y))
                bfs(bfs_write, :) = [temp_x-1, temp_y];
                object_mask(temp_x-1, temp_y, level) = true;
                boolean_mask(temp_x-1, temp_y) = false;
                bfs_write = bfs_write+1;
            end
            if (temp_y>1 && boolean_mask(temp_x, temp_y-1))
                bfs(bfs_write, :) = [temp_x, temp_y-1];
                object_mask(temp_x, temp_y-1, level) = true;
                boolean_mask(temp_x, temp_y-1) = false;
                bfs_write = bfs_write+1;
            end
            if (temp_x<xdim && boolean_mask(temp_x+1, temp_y))
                bfs(bfs_write, :) = [temp_x+1, temp_y];
                object_mask(temp_x+1, temp_y, level) = true;
                boolean_mask(temp_x+1, temp_y) = false;
                bfs_write = bfs_write+1;
            end
            if (temp_y<ydim && boolean_mask(temp_x, temp_y+1))
                bfs(bfs_write, :) = [temp_x, temp_y+1];
                object_mask(temp_x, temp_y+1, level) = true;
                boolean_mask(temp_x, temp_y+1) = false;
                bfs_write = bfs_write+1;
            end
            bfs_read = bfs_read+1;
        end
    end
end


% expand a mask in every direction by the expansion variable
function expanded_mask = expand_mask(mask, expansion)
    [xdim, ydim] = size(mask);
    expanded_mask = mask;
    %expand it horizontally
    for x=1:xdim
        found_row = find(mask(x,:));
        prev = 0;
        first = 0;
        for y_index=1:length(found_row)
            y = found_row(y_index);
            if prev == 0
                first = y;
            elseif y-prev~=1 %we found a gap
                expanded_mask(x, max(1, first-expansion):min(ydim, prev+expansion)) = true;
                first = y;
            end
            prev = y;
        end
        %finish off the last one
        if first~=0
            expanded_mask(x, max(1, first-expansion):min(ydim, prev+expansion)) = true;
        end
    end

    %expand it vertically
    for y=1:ydim
        found_column = find(mask(:,y));
        prev = 0;
        first = 0;
        for x_index=1:length(found_column)
            x = found_column(x_index);
            if prev == 0
                first = x;
            elseif x-prev~=1 %we found a gap
                expanded_mask(max(1, first-expansion):min(xdim, prev+expansion), y) = true;
                first = x;
            end
            prev = x;
        end
        %finish off the last one
        if first~=0
            expanded_mask(max(1, first-expansion):min(xdim, prev+expansion), y) = true;
        end
    end
end


% find how likely it is that two masks are the same object
function rating = sameness(mask, prev_mask, theta, prev_theta)
    % size rating - the larger the better
    size = min(sum(mask(:)), sum(prev_mask(:))) / max(sum(mask(:)), sum(prev_mask(:)));
    
    % direction rating - the larger the better
    dir = abs(theta - prev_theta);
    if dir>pi
        dir = 2*pi - dir;
    end
    dir = (pi - dir) / pi;
    
    % overlap rating - the larger the better
    [prev_mask_offset, ~, ~] = get_overlap(theta, prev_mask, mask, true);
    [mask_offset, ~, ~] = get_overlap(theta, mask, prev_mask, false);
    best_offset1 = sum( sum( mask & prev_mask_offset ) );
    best_offset2 = sum( sum( prev_mask & mask_offset ) );
    overlap = max(best_offset1, best_offset2) / max(sum( mask(:) ), sum( prev_mask(:) ));
    
    % total rating - the larger the better
    rating = dir * overlap;
end


% find the median theta in a mask
function [theta, theta_median, theta_offset] = get_theta_median(object_mask, vx, vy)
    if nargin==3
        theta = atan2(vy, vx);
        mask_theta = theta(object_mask);
    elseif nargin==1
        mask_theta = object_mask;
        theta = mask_theta;
    else
        error('Wrong number of arguments to get_theta_median');
    end
    curr_max = sum(mask_theta<0);
    quadrant = 1;
    if sum(mask_theta>-pi/2 & mask_theta<pi/2)>curr_max
        curr_max = sum(mask_theta>-pi/2 & mask_theta<pi/2);
        quadrant = 2;
    end
    if sum(mask_theta>0)>curr_max
        curr_max = sum(mask_theta>0);
        quadrant = 3;
    end
    if sum(mask_theta<-pi/2 | mask_theta>pi/2)>curr_max
        curr_max = sum(mask_theta<-pi/2 | mask_theta>pi/2);
        quadrant = 4;
    end
    
    if quadrant==1
        theta_offset = pi/2;
    elseif quadrant==2
        theta_offset = 0;
    elseif quadrant==3
        theta_offset = -pi/2;
    elseif quadrant==4
        theta_offset = pi;
    end
    theta = theta + theta_offset;
    theta(theta>pi) = theta(theta>pi)-2*pi;
    theta(theta<-pi) = theta(theta<-pi)+2*pi;
    if nargin==3
        theta_median = median(theta(object_mask));
    elseif nargin==1
        theta_median = median(theta);
    else
        error('Wrong number of arguments to get_theta_median');
    end
end