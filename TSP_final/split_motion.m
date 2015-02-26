load('ready_for_motion_split.mat');
[xdim, ydim] = size(IMG_label);
total_flow = zeros(xdim, ydim);
total_flow(IMG_w+1:end-IMG_w, IMG_w+1:end-IMG_w) = abs(flow.bvx) + abs(flow.bvy);

SP_aggregate_flow = zeros(IMG_K, 1);

for count=1:10
    fprintf('split motion iter %d\n', count);
    for k=1:IMG_K
        SP_aggregate_flow(k) = sum(sum(total_flow(IMG_label==k)));
    end

    mean_aggregate = mean(SP_aggregate_flow);

    [sorted_aggregate, sorted_aggregate_indices] = sort(SP_aggregate_flow);

    k = IMG_K;
    while k>0 && sorted_aggregate(k) > mean_aggregate
        fprintf('splitting %d, k=%d\n', sorted_aggregate_indices(k), k);
        [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = move_split_SP(IMG_label, IMG_SP, IMG_K, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_max_SPs, IMG_data, IMG_boundary_mask, model_order_params, false(size(IMG_SP_old)), IMG_alive_dead_changed, IMG_N, IMG_new_SP, true(size(IMG_SP_changed)), IMG_K, sorted_aggregate_indices(k), true);
        k = k-1;
    end
end

[IMG_K, IMG_label, IMG_SP, SP_changed1, IMG_max_UID, IMG_alive_dead_changed, IMG_SP_old] = localonly_move(IMG_label, IMG_K, IMG_N, IMG_SP_changed, IMG_SP, IMG_T4Table, IMG_boundary_mask, IMG_dummy_log_prob, IMG_new_SP, IMG_SP_old, IMG_data, model_order_params, IMG_new_pos, IMG_new_app, IMG_max_UID, IMG_alive_dead_changed, IMG_max_SPs, 150);

%percentage = floor(IMG_K/4);

% for k_index=1:percentage
%     k = sorted_k(k_index);
%     [IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed] = merge_move_internal(IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed, IMG_SP_old, model_order_params, k, true);
% end
% for k_index=IMG_K-percentage:IMG_K
%     k = sorted_k(k_index);
% end