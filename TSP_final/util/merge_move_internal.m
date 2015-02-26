function [IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed] = merge_move_internal(IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed, IMG_SP_old, model_order_params, k, forced)
    % find all bordering super pixels
    neighbors = U_find_border_SP(IMG_label, IMG_SP, k, false(IMG_K, 1));

    max_E = -inf;
    max_k = -1;

    % loop through all neighbors
    found_neighbors = find(neighbors);
    for merge_k_index=1:length(found_neighbors);
        merge_k = found_neighbors(merge_k_index);
        new_E = U_move_merge_calc_delta(IMG_SP(k), IMG_SP(merge_k), IMG_SP_old(k), IMG_SP_old(merge_k), model_order_params);
        if new_E > max_E || max_k==-1
            max_E = new_E;
            max_k = merge_k;
        end
    end

    % merge if it increases energy or if we're forced to merge
    if max_E>0 || (forced && max_k~=-1)
        % change the labels
        IMG_label(IMG_label==max_k) = k;
        IMG_SP_changed(k) = true;
        IMG_SP_changed(max_k) = true;

        IMG_SP = U_merge_SPs(IMG_SP, IMG_label, k, max_k);

        if IMG_SP_old(max_k)
            IMG_alive_dead_changed = true;
        elseif max_k<numel(IMG_SP)
            IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(numel(IMG_SP)), max_k);
            IMG_SP = U_merge_SPs(IMG_SP, IMG_label, max_k, numel(IMG_SP));
            IMG_SP(numel(IMG_SP)) = [];
            IMG_K = IMG_K-1;
        end
    end
end