root = 'sequences/wheelock_small/'; % folder holding flow_folder and image files
flow_folder = 'flow/'; % name of the folder holding the optical flow files in the root directory
input_video = 'video/wheelock_small.mov'; % video to track motion of
output_video = 'results/wheelock_small.avi'; % name of video to create
img_suffix = 'jpg'; % type of images used

params.frame_inc = 5; % interval of frames to look at
params.expansion = 40; % the amount (radius) to expand each object mask by
params.percentage_tolerance = 8; % tolerance for flow direction in radians
params.min_object_percentage = 1; % the minimum percentage of the image an object can be for it to be taken into consideration
params.cutoff_iters = 2; % the number of times to zero out everything below the mean flow
params.sameness_threshold = 0.4; % the sameness amount necessary for two masks to be considered the same object
params.std_threshold = 0.2; % the threshold for a mask to be considered an object
params.min_frames_per_object = 3; % the minimum number of frames that an object can consist of
params.frame_rate = 10; % default frame rate in case we're not reading in from a video

if exist(input_video, 'file')
    params.frame_rate = write_frames(input_video, root, img_suffix);
end

compute_flows(root, img_suffix, flow_folder, params.frame_inc);

ims = get_masked_images(params, root, flow_folder, output_video);

write_video(ims, output_video, params.frame_rate);
