% =============================================================================
% == localonly_move.m
% == --------------------------------------------------------------------------
% == An interface to perform local TSP moves without flow estimation.
% == See m files for calling convention.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-4-2014
% =============================================================================

function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = localonly_move(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs, its)
    disp('localonly_move');
    for iter=1:its
        if mod(iter,10)==0
            fprintf('localonly iter=%d\n',iter);
        end
        imagesc(IMG_label);
        drawnow;
        disp('drawing');
        [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old, changed] = local_move_internal(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs);
        if ~changed
            break;
        end
    end
end
