function compute_flows(root, img_suffix, flow_folder, frame_inc)
    addpath('util/');
    files = dir([root '*.' img_suffix]);
    root_flows = fullfile(root, flow_folder);
    if (~exist(root_flows,'dir'))
        mkdir(root_flows);
    end

    disp('Precomputing all the optical flows...');
    for f=(1+frame_inc*2):frame_inc:numel(files)
        im1 = imread(fullfile(root,files(f-frame_inc*2).name));
        im2 = imread(fullfile(root,files(f).name));
        outname = fullfile(root_flows,[files(f).name(1:end-4) '_flow.mat']);
        disp([' -> ' outname]);
        compute_of(im1,im2,outname);
    end
    disp(' -> Optical flow calculations done');
end