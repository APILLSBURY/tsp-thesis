% =============================================================================
% == split_move.cpp
% == --------------------------------------------------------------------------
% == A MEX interface to perform split moves on TSPs.
% == See m files for calling convention.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-05-2014
% =============================================================================


function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = split_move(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, IMG_N, IMG_new_SP, its)
    disp('split_move');
    for i=1:its
        % choose a random order of super pixels
        Nsp = numel(IMG_SP);
        perm = randperm(Nsp);

        pre_K = Nsp;
        %split_thres = floor(IMG_N/IMG_K);

        energies = zeros(Nsp, 1);
        mean_area = 0;

        for k=1:pre_K
            SP_size = IMG_SP(k).N;
            if SP_size==0
                model_order = 0;
            else
                size_var = (model_order_params.area - SP_size/2)*SP_size/model_order_params.area_var;
                if IMG_SP_old(k)
                    model_order = model_order_params.is_old_const + size_var;
                else
                    model_order = model_order_params.is_new_const + size_var;
                end
            end
            
            temp_energy = (IMG_SP(k).log_likelihood + model_order) / IMG_SP(k).N;
            energies(k) = temp_energy;
            mean_area = mean_area + IMG_SP(k).N;
        end
        mean_area = mean_area / pre_K;
        threshold = min(energies) + (max(energies)-min(energies))*0.2;

        for kindex=1:length(perm)
            k = perm(kindex);
            if k<=numel(IMG_SP) && (any(IMG_SP(k).neighbors) && (IMG_SP(k).N>mean_area || energies(k) < threshold) && IMG_SP_changed(k))
                [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = move_split_SP(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_SP_old, IMG_alive_dead_changed, IMG_N, IMG_new_SP, true(size(IMG_SP_changed)), Nsp, k);
            end
        end
    end
end