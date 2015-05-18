function [final_masks] = concatenate_masks()
    [masks, theta_medians] = find_masks();
    final_masks = false(size(masks));
    for frame=2:size(masks, 3)-1
        reverse_theta = theta_medians(frame) + pi;
        if reverse_theta>pi
            reverse_theta = reverse_theta - 2*pi;
        end
        before_overlap = get_overlap(reverse_theta, masks(:,:,frame-1), masks(:,:,frame));
        after_overlap = get_overlap(theta_medians(frame+1), masks(:,:,frame+1), masks(:,:,frame));
        
        %final_masks(:,:,frame) = (before_overlap & masks(:,:,frame)) | (before_overlap & after_overlap) | (masks(:,:,frame) & after_overlap);
        final_masks(:,:,frame) = before_overlap | masks(:,:,frame) | after_overlap;
    end
    
    %apply the mask to the image
    im = imread('sequences/every_5/00016.jpg');
    masked_im = im;
    masked_im(:,:,1) = uint8(double(im(:,:,1)) .* double(final_masks(:,:,2)));
    masked_im(:,:,2) = uint8(double(im(:,:,2)) .* double(final_masks(:,:,2)));
    masked_im(:,:,3) = uint8(double(im(:,:,3)) .* double(final_masks(:,:,2)));
    image(masked_im);
end

function [max_offset_mask] = get_overlap(theta_median, moving_mask, stationary_mask)
    [xdim, ydim] = size(moving_mask);
    max_overlap = -1;
    max_offset_mask = false(xdim, ydim);
    rho = 0;
    while rho==0 || sum(offset_mask(:))>0
        [col_offset, row_offset] = pol2cart(theta_median, rho);
        row_offset = round(row_offset);
        col_offset = round(col_offset);
        offset_mask = false(xdim, ydim);
        offset_mask(max(1+row_offset, 1):min(end, end+row_offset), max(1+col_offset, 1):min(end, end+col_offset)) = moving_mask(max(1, 1-row_offset):min(end-row_offset, end), max(1, 1-col_offset):min(end-col_offset, end));
        if mod(rho, 20)==1
            imagesc(offset_mask | stationary_mask);
            drawnow;
        end
        combined_mask = offset_mask & stationary_mask;
        overlap = sum(combined_mask(:));
        if overlap>max_overlap
            max_overlap = overlap;
            max_offset_mask = offset_mask;
        end
        rho = rho + 1;
    end
end

function [masks, theta_medians] = find_masks()
    FIRST_FRAME = 11; % number of the frame to start at
    LAST_FRAME = 21; % number of the frame to end at
    FRAME_INC = 5; % amount to increment the frame numbers by
    EXPANSION = 20; % the amount (radius) to expand each object mask by
    PERCENTAGE_TOLERANCE = 2; % tolerance for flow direction in radians
    MIN_OBJECT_PERCENTAGE = 0.1; % the minimum percentage of the image an object can be for it to be taken into consideration
    CUTOFF_ITERS = 4; % the number of times to zero out everything below the mean flow
    FOLDER = 'every_5';
    
    TOLERANCE = 2*pi*PERCENTAGE_TOLERANCE/100;
    NUM_FRAMES = floor((LAST_FRAME-FIRST_FRAME)/FRAME_INC)+1;
    
    for frame=1:NUM_FRAMES
        frame_num = (frame-1)*FRAME_INC + FIRST_FRAME;
        % flow from past
        load(['sequences/' FOLDER '/TSP_flows/000' num2str(frame_num) '_flow.mat']);
        vx1 = flow.bvx;
        vy1 = flow.bvy;
        [xdim, ydim] = size(vx1);
        additive_flow = abs(vx1)+abs(vy1);

        % create the masks matrix
        if frame==1
            masks = zeros(xdim, ydim, NUM_FRAMES);
            theta_medians = zeros(NUM_FRAMES);
        end
        
        %make a mask of the flow values
        additive_flow = additive_flow ./ max(additive_flow(:));
        for iter=1:CUTOFF_ITERS
            cutoff = mean(mean(additive_flow(additive_flow>0)));
            additive_flow(additive_flow<cutoff) = 0;
        end
        additive_flow(additive_flow>=cutoff) = 1;

        boolean_mask = false(size(additive_flow));
        boolean_mask(additive_flow==1) = true;

        % make a separate mask for each object
        object_mask = false(xdim, ydim, 10);
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
            level = level + 1;
            bfs(bfs_write, :) = [temp_x, temp_y];
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

        for l=1:level
            %expand the mask to make sure we don't miss anything
            expanded_mask = object_mask(:,:,l);
            if length(find(expanded_mask)) < (xdim*ydim)*(MIN_OBJECT_PERCENTAGE/100)
                continue;
            end

            %expand it horizontally
            for x=1:xdim
                found_row = find(object_mask(x,:,l));
                prev = 0;
                first = 0;
                for y_index=1:length(found_row)
                    y = found_row(y_index);
                    if prev == 0
                        first = y;
                    elseif y-prev~=1 %we found a gap
                        expanded_mask(x, max(1, first-EXPANSION):min(ydim, prev+EXPANSION)) = true;
                        first = y;
                    end
                    prev = y;
                end
                %finish off the last one
                if first~=0
                    expanded_mask(x, max(1, first-EXPANSION):min(ydim, prev+EXPANSION)) = true;
                end
            end

            %expand it vertically
            for y=1:ydim
                found_column = find(object_mask(:,y,l));
                prev = 0;
                first = 0;
                for x_index=1:length(found_column)
                    x = found_column(x_index);
                    if prev == 0
                        first = x;
                    elseif x-prev~=1 %we found a gap
                        expanded_mask(max(1, first-EXPANSION):min(xdim, prev+EXPANSION), y) = true;
                        first = x;
                    end
                    prev = x;
                end
                %finish off the last one
                if first~=0
                    expanded_mask(max(1, first-EXPANSION):min(xdim, prev+EXPANSION), y) = true;
                end
            end

            [theta, theta_median, theta_offset] = get_theta_median(object_mask(:,:,l), vx1, vy1);

            % make sure that we're accounting for the values wrapping around
            if theta_median < TOLERANCE-pi
                expanded_mask(theta>theta_median+TOLERANCE & theta<theta_median-TOLERANCE+2*pi) = false;
            elseif theta_median > pi-TOLERANCE
                expanded_mask(theta>theta_median+TOLERANCE-2*pi & theta<theta_median-TOLERANCE) = false;
            else
                expanded_mask(theta>theta_median+TOLERANCE | theta<theta_median-TOLERANCE) = false;
            end

            % check to see if this is the biggest object
            if sum(expanded_mask(:)) > sum(sum(masks(:,:,frame)))
                theta_median = theta_median - theta_offset;
                if theta_median>pi
                    theta_median = theta_median - 2*pi;
                end
                if theta_median<-pi
                    theta_median = theta_median + 2*pi;
                end
                theta_medians(frame) = theta_median;
                masks(:,:,frame) = expanded_mask;
            end
        end
    end
end

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

function write_video(frames, vid_name)
% Create a VideoWriter object to write the video out to a new, different file.
	writerObj = VideoWriter(vid_name);
	open(writerObj);
	
	% Read the frames back in from disk, and convert them to a movie.
	% Preallocate recalledMovie, which will be an array of structures.
	% First get a cell array with all the frames.
	allTheFrames = cell(numberOfFrames,1);
	allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
	% Next get a cell array with all the colormaps.
	allTheColorMaps = cell(numberOfFrames,1);
	allTheColorMaps(:) = {zeros(256, 3)};
	% Now combine these to make the array of structures.
	recalledMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps)
	for frame = 1 : numberOfFrames
		% Construct an output image file name.
		outputBaseFileName = sprintf('Frame %4.4d.png', frame);
		outputFullFileName = fullfile(outputFolder, outputBaseFileName);
		% Read the image in from disk.
		thisFrame = imread(outputFullFileName);
		% Convert the image into a "movie frame" structure.
		recalledMovie(frame) = im2frame(thisFrame);
		% Write this frame out to a new video file.
		writeVideo(writerObj, thisFrame);
	end
	close(writerObj);
end