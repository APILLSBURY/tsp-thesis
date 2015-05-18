addpath('motion_tracker/');
root = 'sequences/tuck/'; % folder holding flow_folder and image files
flow_folder = 'flow/'; % name of the folder holding the optical flow files in the root directory
input_video = 'tuck.avi'; % video to track motion of
output_video = 'results/tuck.avi'; % name of video to create
img_suffix = 'jpg'; % type of images used
frame_inc = 5; % interval of frames to look at
expansion = 40; % the amount (radius) to expand each object mask by
percentage_tolerance = 8; % tolerance for flow direction in radians
min_object_percentage = 1; % the minimum percentage of the image an object can be for it to be taken into consideration
cutoff_iters = 2; % the number of times to zero out everything below the mean flow
test = false; % whether we're testing on a specific frame
sameness_threshold = 0.4;
std_threshold = 0.2;
min_frames_per_object = 3;
frame_rate = 10; % default frame rate in case we're not reading in from a video

if test
    frame_inc = 0;
else
    if exist(input_video, 'file')
        frame_rate = write_frames(input_video, root, img_suffix);
    end
    compute_flows(root, img_suffix, flow_folder, frame_inc);
end

ims = get_masked_images(frame_inc, root, flow_folder, expansion, percentage_tolerance, min_object_percentage, cutoff_iters, test, sameness_threshold, std_threshold, min_frames_per_object);

if test
    if size(ims, 4) > 1
        image(ims(:,:,:,floor(size(ims, 4)/2)));
    else
        image(ims);
    end
else
    write_video(ims, output_video, frame_rate);
end
