function IMG = IMG_prop(img, vx, vy, IMG)
    % delete the SPs that are only in the boundary, i.e are occluded
    unique_vals_image = unique(IMG.label(IMG.boundary_mask));
    unique_vals_boundary = unique(IMG.label(~IMG.boundary_mask));
    boundarySPs = setdiff(unique_vals_boundary, unique_vals_image);
    IMG.label(ismember(IMG.label, boundarySPs)) = 0;


    % Build IMG structure

    % 1. image statistics
    IMG.data = SP_img2data(img, IMG.w);
    IMG.data = rescale_data(IMG.data, IMG.hyper.op_Sigma, IMG.hyper.oa_Sigma);

    % DW: hard for later tracking...
    % delete unused super pixels and change the label as well
    labelmap = zeros(IMG.K, 1);
    k = 1;
    totali = 1;
    mask = IMG.label>0;

    % delete empty SPs, relabel the others to make up for it
    while k<=IMG.K
        if SP_is_empty(IMG, k)
            IMG.SP(k) = [];
            IMG.K = IMG.K - 1;
        else
            labelmap(totali) = k;
            k = k + 1;
        end
        totali = totali + 1;
    end
    IMG.label(mask) = labelmap(IMG.label(mask));

    % get the means of the SPs
    IMG.prev_pos_mean = zeros(IMG.K, 2);
    IMG.prev_app_mean = zeros(IMG.K, 3);
    for k=1:IMG.K
        if ~SP_is_empty(IMG, k)
            IMG.prev_pos_mean(k, :) = NormalD_calc_mean(IMG.SP(k).pos);
            IMG.prev_app_mean(k, :) = NormalD_calc_mean(IMG.SP(k).app);
        end
    end
    meanx = IMG.prev_pos_mean(:, 1);
    meany = IMG.prev_pos_mean(:, 2);

    IMG.prev_label = IMG.label;
    IMG.prev_K = IMG.K;



    % get the previous precision
    mu_p = bsxfun(@times, IMG.prev_pos_mean, sqrt(IMG.hyper.op_Sigma));
    mu_a = bsxfun(@times, IMG.prev_app_mean, sqrt(IMG.hyper.oa_Sigma));
    mu = [mu_p, mu_a];
    [IMG.prev_covariance, IMG.prev_precision] = get_gp_covariance(IMG.label, mu, IMG.cov_var_a, IMG.cov_var_p, IMG.hyper.p_Delta(1));

    vx_extended = holdpad(vx, size(IMG.label,1),size(IMG.label,2));
    vy_extended = holdpad(vy, size(IMG.label,1),size(IMG.label,2));

    IMG.prev_indices = populate_indices(IMG.prev_K, IMG.prev_label);

    sp_v = zeros(2, numel(IMG.SP));
    sp_x = zeros(1, numel(IMG.SP));
    sp_y = zeros(1, numel(IMG.SP));
    for k = 1:IMG.K
        indices_k = IMG.prev_indices(k).all;
        IMG.SP(k).app.theta = IMG.prev_app_mean(k, :);

        vxi = mean(vx_extended(indices_k));
        vyi = mean(vy_extended(indices_k));

        x = round((meanx(k))*sqrt(IMG.hyper.op_Sigma(1)))+1;
        y = round((meany(k))*sqrt(IMG.hyper.op_Sigma(2)))+1;
        xi = max(min(x,IMG.xdim-2*IMG.w),1);
        yi = max(min(y,IMG.ydim-2*IMG.w),1);

        %2.2 Pixel statistics
        IMG.SP(k).N = 0;
        IMG.SP(k).pixels = false(size(IMG.SP(k).pixels));
        IMG.SP(k).borders = false(size(IMG.SP(k).borders));
        IMG.SP(k).neighbors = zeros(size(IMG.SP(k).neighbors));
        IMG.SP_old(k) = true;

        IMG.SP(k).pos.theta = IMG.prev_pos_mean(k, :);

        if (isempty(IMG.SP(k).prev_v) || all(IMG.SP(k).prev_v==0))
            IMG.SP(k).prev_v(:) = [vxi, vyi] ./ sqrt(IMG.hyper.op_Sigma);
        else
            IMG.SP(k).prev_v(:) = IMG.SP(k).v(:);
        end

        IMG.SP(k).v = [vx(xi,yi), vy(xi,yi)] ./ sqrt(IMG.hyper.op_Sigma);

        sp_v(1,k) = vxi;
        sp_v(2,k) = vyi;

        sp_x(k) = vxi + x;
        sp_y(k) = vyi + y;
    end
    meanx = meanx * sqrt(IMG.hyper.op_Sigma(1));
    meany = meany * sqrt(IMG.hyper.op_Sigma(2));


    sp_vx = sp_v(1,:);
    sp_vy = sp_v(2,:);

    IMG.label = SP_prop_init(IMG.K, IMG.label, meanx, meany, sp_vx, sp_vy, IMG.boundary_mask);
    IMG.K = max(max(IMG.label));
    IMG = IMG_populate_SPs(IMG);

    IMG.SP_changed(:) = true;
end