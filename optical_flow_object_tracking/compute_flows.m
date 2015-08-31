function compute_flows(root, img_suffix, flow_folder, frame_inc, little_flow)
    addpath('flow_code_v2/');
    files = dir([root '*.' img_suffix]);
    root_flows = fullfile(root, flow_folder);
    if (~exist(root_flows,'dir'))
        mkdir(root_flows);
    end

    % note: this part could be parrallelized to compute flows faster as
    % they don't have to be computed in order
    disp('Precomputing all the optical flows...');
    if little_flow
        parfor f=2:numel(files)
            im1 = imread(fullfile(root,files(f-frame_inc).name));
            im2 = imread(fullfile(root,files(f).name));
            outname = fullfile(root_flows,[files(f).name(1:end-4) '_flow.mat']);
            disp([' -> ' outname]);
            compute_of(im1,im2,outname);
        end
    else
        for f=(1+frame_inc*2):frame_inc:numel(files)
            im1 = imread(fullfile(root,files(f-frame_inc*2).name));
            im2 = imread(fullfile(root,files(f).name));
            outname = fullfile(root_flows,[files(f).name(1:end-4) '_flow.mat']);
            disp([' -> ' outname]);
            compute_of(im1,im2,outname);
        end
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