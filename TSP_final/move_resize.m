function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = move_resize(IMG_label, IMG_w, flow, IMG_K, IMG_SP, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, IMG_alive_dead_changed, IMG_N, IMG_new_SP, IMG_SP_old, iters)
    [xdim, ydim] = size(IMG_label);
    total_flow = zeros(xdim, ydim);
    total_flow(IMG_w+1:end-IMG_w, IMG_w+1:end-IMG_w) = abs(flow.bvx) + abs(flow.bvy);

    denominator = 2;
    for iter=1:iters
        fprintf('split motion iter %d\n', iter);

        % calculate gravity where mass is the SP_aggregate_flow
        disp('calculating gravity');
        SP_mass = zeros(IMG_K, 1);
        for k=1:IMG_K
            SP_mass(k) = sum(sum(total_flow(IMG_label==k)))/IMG_SP(k).N;
        end

        SP_position = zeros(IMG_K, 2);
        for k=1:IMG_K
            if IMG_SP(k).N ~= 0
                SP_position(k, :) =  NormalD_calc_mean(IMG_SP(k).pos);
            end
        end

        SP_gravity = zeros(IMG_K, 1);
        for k1=1:IMG_K
            if ~(k1 > numel(IMG_SP) || isempty(IMG_SP(k1).N) || IMG_SP(k1).N == 0)
                for k2=1:IMG_K
                    if k1~=k2 && ~(k2 > numel(IMG_SP) || isempty(IMG_SP(k2).N) || IMG_SP(k2).N == 0)
                        distance_sq = (SP_position(k1, 1)-SP_position(k2, 1))^2  + (SP_position(k1, 2)-SP_position(k2, 2))^2;
                        SP_gravity(k1) = SP_gravity(k1) + (SP_mass(k1) * SP_mass(k2) / distance_sq);
                    end
                end
            end
        end
        [sorted_gravity, sorted_gravity_indices] = sort(SP_gravity);
        sorted_gravity_indices(sorted_gravity==0) = [];

        fraction = floor(length(sorted_gravity_indices)/denominator);

        for k_index=length(sorted_gravity_indices)-fraction:length(sorted_gravity_indices)
            k = sorted_gravity_indices(k_index);
            if IMG_SP(k).N ~= 0
                if iter <= iters/2
                    fprintf('splitting %d, k=%d\n', k, k_index);
                    [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = move_split_SP(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, false(size(IMG_SP_old)), IMG_alive_dead_changed, IMG_N, IMG_new_SP, true(size(IMG_SP_changed)), IMG_K, k, true);
                else
                    fprintf('merging %d, k=%d\n', k, k_index);
                    [IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed] = merge_move_internal(IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed, IMG_SP_old, model_order_params, k, true);
                end
            end
        end
        if iter <= iters/2
            denominator = denominator+1;
        else
            denominator = denominator*2;
        end
    end
end

