addpath('util/')
root = 'sequences/wheelock3/'; % folder holding flow_folder and image files
flow_folder = 'flow/'; % name of the folder holding the optical flow files in the root directory
fine_flow_folder = 'fine_flow/';
input_video = 'video/wheelock3.mov'; % video to track motion of
output_video = 'results/wheelock3.avi'; % name of video to create
img_suffix = 'jpg'; % type of images used

params.sameness_threshold = 0.4; % the sameness amount necessary for two masks to be considered the same object
params.time_inc = 0.1; % time interval of frames to look at, in seconds
params.expansion = 0; % the amount (radius) to expand each object mask by
params.percentage_tolerance = 6; % tolerance for flow direction in radians
params.min_object_percentage = 0.1; % the minimum percentage of the image an object can be for it to be taken into consideration
params.cutoff_iters = 2; % the number of times to zero out everything below the mean flow
params.min_frames_per_object = round(2 / params.time_inc); % the minimum number of frames that an object can consist of

frame_rate = 10; % default frame rate in case we're not reading in from a video
create_individual_object_videos = true; % whether to create an individual video for each object or not

% write the frames from the input video to the disc
if exist(input_video, 'file')
    frame_rate = write_frames(input_video, root, img_suffix);
else
    error('Input video does not exist.');
end

params.frame_inc = round(params.time_inc * frame_rate);

% make a new set of params for finding fine masks
fine_params = params; 
fine_params.expansion = 40;
fine_params.min_frames_per_object = 10;
fine_params.sameness_threshold = 0.2;

% compute the optical flows for the frames
compute_flows(root, img_suffix, flow_folder, params.frame_inc, false);

% get the masks for each object
object_masks = get_object_masks(params, root, flow_folder);

% compute the fine optical flows for the frames
compute_flows(root, img_suffix, fine_flow_folder, 1, true);

% get the parts of the masks we missed
object_masks = get_fine_masks(object_masks, fine_params, root, fine_flow_folder);

save('object_masks.mat', 'object_masks');

% create videos from the object masks
create_videos(object_masks, root, params.frame_inc, output_video, frame_rate, create_individual_object_videos);
