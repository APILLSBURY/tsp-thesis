root = 'sequences/hop/'; % folder holding flow_folder and image files
flow_folder = 'flow/'; % name of the folder holding the optical flow files in the root directory
input_video = 'video/hop.mp4'; % video to track motion of
output_video = 'results/hop.avi'; % name of video to create
img_suffix = 'jpg'; % type of images used

params.sameness_threshold = 0.4; % the sameness amount necessary for two masks to be considered the same object
params.time_inc = 0.5; % time interval of frames to look at, in seconds
params.expansion = 20; % the amount (radius) to expand each object mask by
params.percentage_tolerance = 6; % tolerance for flow direction in radians
params.min_object_percentage = 0.1; % the minimum percentage of the image an object can be for it to be taken into consideration
params.cutoff_iters = 2; % the number of times to zero out everything below the mean flow
params.std_threshold = 0; % the threshold for a mask to be considered an object
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

% compute the optical flows for the frames
compute_flows(root, img_suffix, flow_folder, params.frame_inc, false);

% get the masks for each object
object_masks = get_object_masks(params, root, flow_folder);

% create videos from the object masks
create_videos(object_masks, root, params.frame_inc, output_video, frame_rate, create_individual_object_videos);
