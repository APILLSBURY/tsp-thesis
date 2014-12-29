% =============================================================================
% == local_move.m
% == --------------------------------------------------------------------------
% == An interface to perform local TSP moves with flow estimation.
% == See m files for calling convention.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-4-2014
% =============================================================================
function IMG = local_move(IMG, its)
    addpath('mex/class_utils/');
    for i=1:its
        IMG = Flow_QP2(IMG);
        IMG = local_move_internal(IMG);
        if ~IMG.changed
            break;
        end
    end
end


function IMG = Flow_QP2(IMG)
    %given label z, globally update flow variable of SP using GP approximation of Optical Flow like L2 penalty

    % get the alive stuff
    num_alive = 1;
    alive2all = zeros(IMG.prev_K, 1);
    for k=1:IMG.prev_K
        if ~SP_is_empty(IMG, k)
            alive2all(num_alive) = k;
            num_alive = num_alive + 1;
        end
    end

    % for computational efficiency
    % 1. Covariance matrix
    if IMG.alive_dead_changed || SxySyy==null ???
        Syy = zeros(num_alive, num_alive);
        Sxy = zeros(IMG.prev_K, num_alive);
        
        %obs_uv(:,1) is obs_u, obs_uv(:,2) is obs_v
        obs_uv = zeros(num_alive,2);

        [Sxy, Syy, SxySyy] = find_SxySyy(all2alive, alive2all, all2dead, dead2all);
        IMG.alive_dead_changed = false;
    end

    % populate the observation
    for k1=1:num_alive
        all_k1 = alive2all(k1);
        obs_uv(k1,:) =IMG.SP(all_k1).pos.mean - IMG.prev_pos_mean(all_k1) - IMG.SP(all_k1).prev_v;
    end

    % uv_gsl(1) = ugsl, uv_gsl(2) = vgsl
    uv_gsl = Sxy * Syy * obs_uv;

    % copy over the new parameters
    for i=1:IMG.prev_k
        if (~isempty(IMG.SP(i).N))
            % ith element of ugsl
            IMG.SP(i).pos.offset = uv_gsl(i,:) + IMG.SP(i).prev_v;
        end
    end
    IMG.SP_changed(1:IMG.prev_k) = true;
end