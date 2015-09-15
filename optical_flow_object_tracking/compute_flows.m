function compute_flows(root, img_suffix, flow_folder, frame_inc, fine_flow)
    addpath('flow_code_v2/');
    files = dir([root '*.' img_suffix]);
    root_flows = fullfile(root, flow_folder);
    if (~exist(root_flows,'dir'))
        mkdir(root_flows);
    end

    if fine_flow
        frame_diff = frame_inc;
    else
        frame_diff = frame_inc*2;
    end
    start_index = 1 + frame_diff;
    
    disp('Precomputing all the optical flows...');
    parfor ff = 0:floor((numel(files)-start_index)/frame_inc)
        f = ff*frame_inc + start_index;
        im1 = imread(fullfile(root,files(f-frame_diff).name));
        im2 = imread(fullfile(root,files(f).name));
        outname = fullfile(root_flows,[files(f).name(1:end-4) '_flow.mat']);
        disp([' -> ' outname]);
        compute_of(im1,im2,outname);
    end
    
    disp(' -> Optical flow calculations done');
end

function [] = compute_of(oim1, oim2, outname)
    if ~exist(outname, 'file')
        [~, bv] = evalc('estimate_flow_interface(oim2, oim1, ''classic+nl-fast'');');
        [~, fv] = evalc('estimate_flow_interface(oim1, oim2, ''classic+nl-fast'');');
        flow.bvx = bv(:, :, 1);
        flow.bvy = bv(:, :, 2);
        flow.fvx = fv(:, :, 1);
        flow.fvy = fv(:, :, 2);
        save(outname, 'flow');
    end
end