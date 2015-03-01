function [sp_labels] = TSP(K, root, files, dispOn, frames)
    %TSP Temporal Superpixel Segmentation.
    %   SP_LABELS = TSP(K, ROOT, FILES) returns the label matrix in time and
    %   space for the video volume in UINT32. K is the (approximate) number of
    %   superpixels per frame. ROOT is the directory to the frames. FILES is a
    %   list of the frame images, typically obtained using
    %      FILES = dir([ROOT '*.jpg']);
    %
    %   SP_LABELS = TSP(K, ROOT, FILES, DISPON) supplies an additional flag to
    %   display the progress of the algorithm while processing. If omitted or
    %   empty, DISPON defaults to true.
    %
    %   SP_LABELS = TSP(K, ROOT, FILES, DISPON, FRAMES) supplies an additional
    %   variable that indicates which frames to process. FRAMES should be in
    %   the format of STARTFRAME:ENDFRAME. If omitted or empty, FRAMES defaults
    %   to 1:NUMFRAMES.

    %   Notes: This version of the code does not reestimate the flow between
    %   frames. As noted in the paper, the flow estimation does not do much. If
    %   you desire to reestimate the flow, set params.reestimateFlow to be
    %   true.
    %
    %   All work using this code should cite:
    %   J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
    %      Temporal Superpixels. CVPR 2013.
    %
    %   Written by Jason Chang and Donglai Wei 2013/06/20

    % add the necessary paths
    addpath('gui/');
    addpath('util/');
    addpath('IMG_utils/');


    IMG_params.cov_var_p = 1000;
    IMG_params.cov_var_a = 100;
    IMG_params.area_var = 400;
    IMG_params.alpha = -15;
    IMG_params.beta = -10;
    IMG_params.deltap_scale = 1e-3;
    IMG_params.deltaa_scale = 100;
    IMG_params.K = K;
    IMG_params.Kpercent = 0.5;
    IMG_params.reestimateFlow = false;

    if (~exist('dispOn','var') || isempty(dispOn))
        dispOn = true;
    end

    root_flows = fullfile(root,'TSP_flows/');
    if (~exist(root_flows,'dir'))
        mkdir(root_flows);
    end

    disp('Precomputing all the optical flows...');
    for f=2:numel(files)
        im1 = imread(fullfile(root,files(f-1).name));
        im2 = imread(fullfile(root,files(f).name));
        outname = fullfile(root_flows,[files(f).name(1:end-4) '_flow.mat']);
        disp([' -> ' outname]);
        compute_of(im1,im2,outname);
    end
    disp(' -> Optical flow calculations done');    


    flow_files = dir([root_flows '*_flow.mat']);

    if (~exist('frames','var') || isempty(frames))
        frames = 1:numel(files);
    else
        frames(frames>numel(files)) = [];
    end
    oim = imread([root files(1).name]);
    sp_labels = zeros(size(oim,1), size(oim,2), numel(frames));
    frame_it = 0;

    disp('Starting Segmentation');
    for f=frames
        disp([' -> Frame '  num2str(f) ' / ' num2str(numel(frames))]);

        frame_it = frame_it + 1;
        oim1 = imread([root files(f).name]);

        if (frame_it==1)
            
            % INITIALIZE IMAGE
            disp('initializing IMG');
     
            IMG_cov_var_a = IMG_params.cov_var_a;
            IMG_cov_var_p = IMG_params.cov_var_p;
            IMG_alive_dead_changed = true;

            % Build IMG structure

            % 1. image statistics
            IMG_oxdim = size(oim1,1);
            IMG_oydim = size(oim1,2);

            N = IMG_oxdim*IMG_oydim;
            IMG_area = N/IMG_params.K;
            IMG_area_var = IMG_params.area_var;


            IMG_w = round(2*sqrt(IMG_area/pi)); %used to be *2

            IMG_xdim = IMG_oxdim + 2*IMG_w;
            IMG_ydim = IMG_oydim + 2*IMG_w;

            IMG_boundary_mask = false(IMG_xdim, IMG_ydim);
            IMG_boundary_mask(IMG_w+1:end-IMG_w, IMG_w+1:end-IMG_w) = true;
            IMG_N = IMG_xdim * IMG_ydim;


            IMG_log_alpha = IMG_params.alpha * IMG_area;
            IMG_log_beta = IMG_params.beta * IMG_area;
            Sigma = IMG_area^2 / (-12.5123*IMG_log_alpha);


            IMG_hyper.p_Sigma = [Sigma Sigma];
            IMG_hyper.p_Delta = [Sigma*2 Sigma*2]*IMG_params.deltap_scale;
            IMG_hyper.a_Sigma = [Sigma*2 Sigma Sigma]*IMG_params.K/100;
            IMG_hyper.a_Delta = [Sigma*20 Sigma*10 Sigma*10]/IMG_params.deltaa_scale;

            r = sqrt(IMG_area / pi);
            IMG_dummy_log_prob = (-0.5 * r^2/IMG_hyper.p_Sigma(1)) - log(2*pi .* IMG_hyper.p_Sigma(1));

            IMG_hyper.p_theta = [0 0];
            IMG_hyper.a_theta = [0 0 0];
            IMG_hyper.op_Sigma = IMG_hyper.p_Sigma;
            IMG_hyper.oa_Sigma = IMG_hyper.a_Sigma;

            IMG_data = SP_img2data(oim1, IMG_w);
            IMG_data = rescale_data(IMG_data, IMG_hyper.op_Sigma, IMG_hyper.oa_Sigma);

            IMG_hyper.p_theta = IMG_hyper.p_theta ./ sqrt(IMG_hyper.p_Sigma);
            IMG_hyper.p_Delta = IMG_hyper.p_Delta ./ (IMG_hyper.p_Sigma);
            IMG_hyper.p_Sigma(:) = 1;

            IMG_hyper.a_theta = IMG_hyper.a_theta ./ sqrt(IMG_hyper.a_Sigma);
            IMG_hyper.a_Delta = IMG_hyper.a_Delta ./ (IMG_hyper.a_Sigma);
            IMG_hyper.a_Sigma(:) = 1;

            % 2. SuperPixel statistics
            IMG_K = round(IMG_params.K*IMG_params.Kpercent);
            IMG_label = random_init(IMG_xdim, IMG_ydim, IMG_w, IMG_K);
            IMG_label(~IMG_boundary_mask) = 0;
            IMG_max_SPs = max(IMG_K + IMG_area, IMG_K * 20);
            
            % 3. Topology Table Look Up
            load topology_tables;
            IMG_T4Table = tc.T4Table;

            IMG_max_UID = 1;

            IMG_SP_changed = true(1, IMG_max_SPs);

            IMG_new_pos = new_NormalD(2, IMG_hyper.p_theta, IMG_hyper.p_Delta, true);
            IMG_new_app = new_NormalD(3, IMG_hyper.a_theta, IMG_hyper.a_Delta, true);
            IMG_new_SP = new_SP(IMG_new_pos, IMG_new_app, 0, [0, 0], IMG_N, IMG_max_SPs);
            IMG_SP_old = false(1,IMG_max_SPs);
            IMG_log_area_var = log(IMG_area_var);

            IMG_K = max(max(IMG_label));
            for i=1:IMG_K
                IMG_SP(i) = new_SP(IMG_new_pos, IMG_new_app, IMG_max_UID, [0, 0], IMG_N, IMG_max_SPs);
                IMG_max_UID = IMG_max_UID+1;
            end

            %CALCULATE VALUES FOR CALCULATING MODEL ORDER
            %calculate model order with is_old/new_const + (area - size/2)*size/area_var
            model_order_params.is_old_const = IMG_log_beta  - 0.5 * (1.837877066409 + IMG_log_area_var + (IMG_area^2)/IMG_area_var);
            model_order_params.is_new_const = IMG_log_alpha - 0.5 * (1.837877066409 + IMG_log_area_var + (IMG_area^2)/IMG_area_var);
            model_order_params.area = IMG_area;
            model_order_params.area_var = IMG_area_var;
            
            %POPULATE SUPERPIXELS
            % populate the linked lists and pointers
            for x=1:IMG_xdim
                for y=1:IMG_ydim
                    curLabel = IMG_label(x, y);
                    if (curLabel>0)
                        index = get_index_from_x_and_y(x, y, IMG_xdim);
                        IMG_SP(curLabel) = SP_add_pixel_init(IMG_SP(curLabel), IMG_data, index, U_check_border_pix(IMG_label, index), IMG_boundary_mask(x, y));
                        IMG_SP(curLabel) = SP_update_neighbors_add_self(IMG_SP(curLabel), IMG_label, index);
                    end
                end
            end

            for k=1:IMG_K
                IMG_SP(k) = SP_calculate_log_probs(IMG_SP(k));
            end
            
            disp('done IMG_init');
        else
            % load the optical flow
            load([root_flows flow_files(f-1).name]);

            vx = -flow.bvx;
            vy = -flow.bvy;
            disp('propagate SPs');
            % IMAGE PROPAGATION
            
            % delete the SPs that are only in the boundary, i.e are occluded
            unique_vals_image = unique(IMG_label(IMG_boundary_mask));
            unique_vals_boundary = unique(IMG_label(~IMG_boundary_mask));
            boundarySPs = setdiff(unique_vals_boundary, unique_vals_image);
            IMG_label(ismember(IMG_label, boundarySPs)) = 0;


            % Build IMG structure

            % 1. image statistics
            IMG_data = SP_img2data(oim1, IMG_w);
            IMG_data = rescale_data(IMG_data, IMG_hyper.op_Sigma, IMG_hyper.oa_Sigma);

            % DW: hard for later tracking...
            % delete unused super pixels and change the label as well
            labelmap = zeros(IMG_K, 1);
            k = 1;
            totali = 1;
            mask = IMG_label>0;

            % delete empty SPs, relabel the others to make up for it
            while k<=IMG_K
                if (k > numel(IMG_SP) || isempty(IMG_SP(k).N) || IMG_SP(k).N == 0)
                    IMG_SP(k) = [];
                    IMG_K = IMG_K - 1;
                else
                    labelmap(totali) = k;
                    k = k + 1;
                end
                totali = totali + 1;
            end
            IMG_label(mask) = labelmap(IMG_label(mask));

            % get the means of the SPs
            IMG_prev_pos_mean = zeros(IMG_K, 2);
            IMG_prev_app_mean = zeros(IMG_K, 3);
            for k=1:IMG_K
                if ~(k > numel(IMG_SP) || isempty(IMG_SP(k).N) || IMG_SP(k).N == 0)
                    IMG_prev_pos_mean(k, :) = NormalD_calc_mean(IMG_SP(k).pos);
                    IMG_prev_app_mean(k, :) = NormalD_calc_mean(IMG_SP(k).app);
                end
            end
            meanx = IMG_prev_pos_mean(:, 1);
            meany = IMG_prev_pos_mean(:, 2);

            IMG_prev_label = IMG_label;
            IMG_prev_K = IMG_K;

            % get the previous precision
            mu_p = bsxfun(@times, IMG_prev_pos_mean, sqrt(IMG_hyper.op_Sigma));
            mu_a = bsxfun(@times, IMG_prev_app_mean, sqrt(IMG_hyper.oa_Sigma));
            mu = [mu_p, mu_a];
            [IMG_prev_covariance, IMG_prev_precision] = get_gp_covariance(mu, IMG_cov_var_a, IMG_cov_var_p, IMG_hyper.p_Delta(1));

            vx_extended = holdpad(vx, size(IMG_label,1),size(IMG_label,2));
            vy_extended = holdpad(vy, size(IMG_label,1),size(IMG_label,2));

            IMG_prev_indices = populate_indices(IMG_prev_K, IMG_prev_label);

            sp_v = zeros(2, numel(IMG_SP));
            sp_x = zeros(1, numel(IMG_SP));
            sp_y = zeros(1, numel(IMG_SP));
            for k = 1:IMG_K
                indices_k = IMG_prev_indices(k).all;
                IMG_SP(k).app.theta = IMG_prev_app_mean(k, :);

                vxi = mean(vx_extended(indices_k));
                vyi = mean(vy_extended(indices_k));

                x = round((meanx(k))*sqrt(IMG_hyper.op_Sigma(1)))+1;
                y = round((meany(k))*sqrt(IMG_hyper.op_Sigma(2)))+1;
                xi = max(min(x,IMG_xdim-2*IMG_w),1);
                yi = max(min(y,IMG_ydim-2*IMG_w),1);

                %2.2 Pixel statistics
                IMG_SP(k).N = 0;
                IMG_SP(k).pixels = false(size(IMG_SP(k).pixels));
                IMG_SP(k).borders = false(size(IMG_SP(k).borders));
                IMG_SP(k).neighbors = zeros(size(IMG_SP(k).neighbors));
                IMG_SP_old(k) = true;

                IMG_SP(k).pos.theta = IMG_prev_pos_mean(k, :);

                if (isempty(IMG_SP(k).prev_v) || all(IMG_SP(k).prev_v==0))
                    IMG_SP(k).prev_v(:) = [vxi, vyi] ./ sqrt(IMG_hyper.op_Sigma);
                else
                    IMG_SP(k).prev_v(:) = IMG_SP(k).v(:);
                end

                IMG_SP(k).v = [vx(xi,yi), vy(xi,yi)] ./ sqrt(IMG_hyper.op_Sigma);

                sp_v(1,k) = vxi;
                sp_v(2,k) = vyi;

                sp_x(k) = vxi + x;
                sp_y(k) = vyi + y;
            end
            meanx = meanx * sqrt(IMG_hyper.op_Sigma(1));
            meany = meany * sqrt(IMG_hyper.op_Sigma(2));


            sp_vx = sp_v(1,:);
            sp_vy = sp_v(2,:);

            IMG_label = SP_prop_init(IMG_K, IMG_label, meanx, meany, sp_vx, sp_vy, IMG_boundary_mask);
            IMG_K = max(max(IMG_label));
            
            % populate the linked lists and pointers
            for x=1:IMG_xdim
                for y=1:IMG_ydim
                    curLabel = IMG_label(x, y);
                    if (curLabel>0)
                        index = get_index_from_x_and_y(x, y, IMG_xdim);
                        IMG_SP(curLabel) = SP_add_pixel_init(IMG_SP(curLabel), IMG_data, index, U_check_border_pix(IMG_label, index), IMG_boundary_mask(x, y));
                        IMG_SP(curLabel) = SP_update_neighbors_add_self(IMG_SP(curLabel), IMG_label, index);
                    end
                end
            end

            for k=1:IMG_K
                IMG_SP(k) = SP_calculate_log_probs(IMG_SP(k));
            end

            IMG_SP_changed(:) = true;
        end

        oim = oim1;

        E = [];
        it = 1;
        IMG_alive_dead_changed = true;
        IMG_Sxy = [];
        IMG_Syy = [];
        converged = false;
        
        while ~converged && it<=5 && frame_it==1
            fprintf('initial while loop: it=%d\n', it);

            oldK = IMG_K;
            IMG_SP_changed(:) = true;

            if (dispOn)
                display_img(IMG_w, IMG_label, it, oim);
            end
            
            [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = split_move(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, IMG_N, IMG_new_SP, 1);
            E(end+1) = U_calc_energy(IMG_N, IMG_SP, IMG_SP_old, model_order_params, IMG_new_SP);
            converged = IMG_K - oldK < 2;

            if (dispOn)
                display_img(IMG_w, IMG_label, it, oim);
            end
            it = it + 1;
        end
        converged = false;
        if (frame_it>1)
            disp('merging, splitting, switching, resizing, localonlying');
            IMG_SP_changed(:) = true;
            [IMG_K, IMG_label, IMG_SP, ~, ~, IMG_SP_old] = merge_move(IMG_label, IMG_SP, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, model_order_params, IMG_K, 1);
            [IMG_K, IMG_label, IMG_SP, ~, IMG_max_UID, ~, IMG_SP_old] = split_move(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, IMG_N, IMG_new_SP, 10);
            [IMG_K, IMG_label, IMG_SP, ~, IMG_max_UID, ~, IMG_SP_old] = switch_move(IMG_label, IMG_SP, IMG_K, IMG_N, IMG_SP_old, IMG_SP_changed, IMG_max_UID, IMG_max_SPs, IMG_alive_dead_changed, IMG_new_SP, model_order_params, IMG_new_pos, IMG_new_app);
            if frame_it==2
                save('pre_resize.mat', IMG_label);
                [IMG_K, IMG_label, IMG_SP, ~, IMG_max_UID, ~, IMG_SP_old] = move_resize(IMG_label, IMG_w, flow, IMG_K, IMG_SP, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_alive_dead_changed, IMG_N, IMG_new_SP, IMG_SP_old, IMG_SP_changed, 10);
            end
            [IMG_K, IMG_label, IMG_SP, ~, IMG_max_UID, ~, IMG_SP_old] = localonly_move(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs, 10);
        end
        IMG_SP_changed(:) = true;
        IMG_alive_dead_changed = true;
        
        it = 0;
        while (~converged && it<20)
            it = it + 1;
            fprintf('\nBig while loop: it=%d\n', it);
            move_times = zeros(1,5);

            if (~IMG_params.reestimateFlow)
                save('pre_localonly.mat');
                disp('localonly_moving');
                tic;[IMG_K, IMG_label, IMG_SP, SP_changed1, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = localonly_move(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs, 150);move_times(2)=toc;
                SP_changed0 = SP_changed1;
            else
                disp('local and localonly');
                tic;[IMG_K, IMG_label, IMG_SP, SP_changed0, IMG_max_UID, IMG_alive_dead_changed, IMG_Sxy, IMG_Syy, IMG_SP_old] = local_move(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_prev_pos_mean, IMG_prev_K, IMG_prev_precision, IMG_prev_covariance, IMG_Sxy, IMG_Syy, IMG_max_SPs, 50);move_times(1)=toc;
                tic;[IMG_K, IMG_label, IMG_SP, SP_changed1, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = localonly_move(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs, 5);move_times(2)=toc;
            end
            if (frame_it>1 && it<5)
                disp('merge, split, switch');
                tic;[IMG_K, IMG_label, IMG_SP, SP_changed2, IMG_alive_dead_changed, IMG_SP_old] = merge_move(IMG_label, IMG_SP, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, model_order_params, IMG_K, 1);move_times(3)=toc;
                tic;[IMG_K, IMG_label, IMG_SP, SP_changed3, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = split_move(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, IMG_N, IMG_new_SP, 1);move_times(4)=toc;
                tic;[IMG_K, IMG_label, IMG_SP, SP_changed4, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = switch_move(IMG_label, IMG_SP, IMG_K, IMG_N, IMG_SP_old, IMG_SP_changed, IMG_max_UID, IMG_max_SPs, IMG_alive_dead_changed, IMG_new_SP, model_order_params, IMG_new_pos, IMG_new_app);move_times(5)=toc;
                IMG_SP_changed = SP_changed0 | SP_changed1 | SP_changed2 | SP_changed3 | SP_changed4;
            else
                IMG_SP_changed = SP_changed0 | SP_changed1;
            end
            
            if ~(frame_it>1 && it<5)
                IMG_SP_changed = SP_changed0 | SP_changed1;
            end

            E(end+1) = U_calc_energy(IMG_N, IMG_SP, IMG_SP_old, model_order_params, IMG_new_SP);
            converged = ~any(~arrayfun(@(x)(isempty(x{1})), {IMG_SP(:).N}) & IMG_SP_changed(1:IMG_K));

            if (dispOn)
                display_img(IMG_w, IMG_label, it, oim);
            end
        end

        SP_UID = {IMG_SP(:).UID};
        found_mask = find(arrayfun(@(x)(isempty(x{1})), SP_UID));
        
        for m = 1:length(found_mask)
            SP_UID{found_mask(m)} = -1;
        end
        sp_labels(:,:,frame_it) = reshape([SP_UID{IMG_label(IMG_w+1:end-IMG_w,IMG_w+1:end-IMG_w)}], size(oim,1), size(oim,2));
    end
end

function display_img(IMG_w, IMG_label, it, oim)
    sfigure(1);
    subplot(1,1,1);
    imagesc(IMG_label);
    title([num2str(it) ' - ' num2str(numel(unique(IMG_label))-1)]);

    sfigure(2);
    subplot(1,1,1);
    im = zeros(size(oim,1)+2*IMG_w, size(oim,2)+2*IMG_w, 3);
    im(IMG_w+1:end-IMG_w, IMG_w+1:end-IMG_w, :) = double(oim)/255;
    borders = is_border_vals(IMG_label);
    im = setPixelColors(im, find(borders), [1 0 0]);
    image(im,'parent',gca);
    drawnow;
end


