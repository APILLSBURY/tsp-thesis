function [] = create_videos(object_masks, root, frame_inc, output_video, frame_rate, create_individual_object_videos)
    % create a video for each object
    if create_individual_object_videos
        for o=1:size(object_masks, 4)
            ims = mask_images(object_masks(:,:,:,o), root, frame_inc);
            split_str = strsplit(output_video, '.');
            output = split_str{1};
            for s=2:(length(split_str)-1)
                output = strcat(output, '.', split_str{s});
            end
            output = strcat(output, '_obj', num2str(o), '.', split_str{end});
            write_video(ims, output, frame_rate);
        end
    end

    %create combined video for all the objects
    combined_final_masks = false(size(object_masks(:,:,:,1)));
    for o=1:size(object_masks, 4)
        combined_final_masks = combined_final_masks | object_masks(:,:,:,o);
    end
    ims = mask_images(combined_final_masks, root, frame_inc);
    
    write_video(ims, output_video, frame_rate);
end


% apply masks to the images
% this puts the masks over the existing images
function ims = mask_images(final_masks, root, frame_inc)
    files = dir([root '*.jpg']);
    ims = uint8(zeros(size(final_masks, 1), size(final_masks, 2), 3, size(final_masks, 3)));
    for f=(2+frame_inc):(size(ims,4)+frame_inc+1)
        findex = f-frame_inc-1;
        ims(:,:,:,findex) = imread(fullfile(root,files(f).name));
        mask = double(final_masks(:,:,findex));
        mask(mask==0) = 0.2;
        
        %apply the mask to the image
        ims(:,:,1,findex) = uint8(double(ims(:,:,1,findex)) .* mask);
        ims(:,:,2,findex) = uint8(double(ims(:,:,2,findex)) .* mask);
        ims(:,:,3,findex) = uint8(double(ims(:,:,3,findex)) .* mask);
    end
end


%this makes black and white images with just the masks
% function ims = mask_images(final_masks, root, frame_inc)
%     ims = uint8(zeros(size(final_masks, 1), size(final_masks, 2), 3, size(final_masks, 3)));
%     for f=(2+frame_inc):(size(ims,4)+frame_inc+1)
%         findex = f-frame_inc-1;
%         mask = double(final_masks(:,:,findex));
%         mask(mask==1) = 255;
%         
%         %apply the mask to the image
%         ims(:,:,1,findex) = uint8(mask);
%         ims(:,:,2,findex) = uint8(mask);
%         ims(:,:,3,findex) = uint8(mask);
%     end
% end

% write a video from a collection of images
function write_video(ims, vid_name, frame_rate)
    disp(['Writing video ' vid_name]);
    % Create a VideoWriter object to write the video out to a new, different file.
	writerObj = VideoWriter(vid_name);
    writerObj.FrameRate = frame_rate;
	open(writerObj);
	for frame = 1 : size(ims, 4)
		writeVideo(writerObj, ims(:,:,:,frame));
	end
	close(writerObj);
end