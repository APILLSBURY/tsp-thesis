function IMG = IMG_init(img, params)

IMG.cov_var_a = params.cov_var_a;
IMG.cov_var_p = params.cov_var_p;
IMG.alive_dead_changed = true;

% Build IMG structure

% 1. image statistics
IMG.oxdim = size(img,1);
IMG.oydim = size(img,2);

N = IMG.oxdim*IMG.oydim;
IMG.area = N/params.K;
IMG.area_var = params.area_var;


IMG.w = round(2*sqrt(IMG.area/pi)); %used to be *2

IMG.xdim = IMG.oxdim + 2*IMG.w;
IMG.ydim = IMG.oydim + 2*IMG.w;

IMG.boundary_mask = false(IMG.xdim, IMG.ydim);
IMG.boundary_mask(IMG.w+1:end-IMG.w, IMG.w+1:end-IMG.w) = true;
IMG.N = IMG.xdim * IMG.ydim;


IMG.log_alpha = params.alpha * IMG.area;
IMG.log_beta = params.beta * IMG.area;
Sigma = IMG.area^2 / (-12.5123*IMG.log_alpha);


IMG.hyper.p_Sigma = [Sigma Sigma];
IMG.hyper.p_Delta = [Sigma*2 Sigma*2]*params.deltap_scale;
IMG.hyper.a_Sigma = [Sigma*2 Sigma Sigma]*params.K/100;
IMG.hyper.a_Delta = [Sigma*20 Sigma*10 Sigma*10]/params.deltaa_scale;

r = sqrt(IMG.area / pi);
IMG.dummy_log_prob = (-0.5 * r^2/IMG.hyper.p_Sigma(1)) - log(2*pi .* IMG.hyper.p_Sigma(1));

IMG.hyper.p_theta = [0 0];
IMG.hyper.a_theta = [0 0 0];
IMG.hyper.op_Sigma = IMG.hyper.p_Sigma;
IMG.hyper.oa_Sigma = IMG.hyper.a_Sigma;

IMG.data = SP_img2data(img, IMG.w);
IMG.data = rescale_data(IMG.data, IMG.hyper.op_Sigma, IMG.hyper.oa_Sigma);

IMG.hyper.p_theta = IMG.hyper.p_theta ./ sqrt(IMG.hyper.p_Sigma);
IMG.hyper.p_Delta = IMG.hyper.p_Delta ./ (IMG.hyper.p_Sigma);
IMG.hyper.p_Sigma(:) = 1;

IMG.hyper.a_theta = IMG.hyper.a_theta ./ sqrt(IMG.hyper.a_Sigma);
IMG.hyper.a_Delta = IMG.hyper.a_Delta ./ (IMG.hyper.a_Sigma);
IMG.hyper.a_Sigma(:) = 1;

% 2. SuperPixel statistics
IMG.K = round(params.K*params.Kpercent);
IMG.label = random_init(IMG.xdim, IMG.ydim, IMG.w, IMG.K);
IMG.label(~IMG.boundary_mask) = 0;
IMG.max_SPs = max(IMG.K + IMG.area, IMG.K * 20);

% 3. Topology Table Look Up
load topology_tables;
IMG.T4Table = tc.T4Table;

IMG.max_UID = 1;

IMG.SP_changed = true(1, IMG.max_SPs);

IMG.new_pos = new_NormalD(2, IMG.hyper.p_theta, IMG.hyper.p_Delta, true);
IMG.new_app = new_NormalD(3, IMG.hyper.a_theta, IMG.hyper.a_Delta, true);
IMG.SP_old = false(1,IMG.max_SPs);
IMG.log_area_var = log(IMG.area_var);

disp('initializing IMG');
IMG.K = max(max(IMG.label));
for i=1:IMG.K
    IMG.SP(i) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, [0, 0], IMG.N, IMG.max_SPs);
    IMG.max_UID = IMG.max_UID+1;
end

IMG = IMG_populate_SPs(IMG);
