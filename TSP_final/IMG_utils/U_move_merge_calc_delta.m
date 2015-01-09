% --------------------------------------------------------------------------
% -- move_merge_calc_delta
% --   calculates the probability of merging the superpixel at k with the
% --   superpixel at merge_k.
% --
% --   parameters:
% --     - k : the original superpixel
% --     - merge_k : the superpixel to merge with
% --------------------------------------------------------------------------

function prob = U_move_merge_calc_delta(IMG_SP_k, IMG_SP_merge_k, k_is_old, merge_k_is_old, model_order_params)
    prob = SP_log_likelihood_test_merge1(IMG_SP_k, IMG_SP_merge_k) + IMG_SP_merge_k.log_likelihood_empty;
    prob = prob - IMG_SP_k.log_likelihood - IMG_SP_merge_k.log_likelihood;
    
%     prob = prob + U_calc_model_order(IMG, IMG_SP_k.N + IMG_SP_merge_k.N, IMG_SP_old_k);
%     prob = prob - U_calc_model_order(IMG, IMG_SP_k.N, IMG_SP_old_k);
%     prob = prob - U_calc_model_order(IMG, IMG_SP_merge_k.N, merge_k_is_old);
    
    k_size = IMG_SP_k.N;
    merge_k_size = IMG_SP_merge_k.N;
    if k_size>0 && merge_k_size>0
        k_and_merge_k_size_var = (model_order_params.area - (k_size+merge_k_size)/2)*(k_size+merge_k_size)/model_order_params.area_var;
        merge_k_size_var = (model_order_params.area - merge_k_size/2)*merge_k_size/model_order_params.area_var;
        k_size_var = (model_order_params.area - k_size/2)*k_size/model_order_params.area_var;
        prob = prob + k_and_merge_k_size_var - k_size_var - merge_k_size_var;
        if merge_k_is_old
            prob = prob + model_order_params.is_old_const;
        else
            prob = prob + model_order_params.is_old_const;
        end
    elseif merge_k_size>0 %if k_size == 0
        if k_is_old~=merge_k_is_old % if they're the same, the constants cancel out when we add one and subtract the other
            if k_is_old % k is old, merge_k is new
                prob = prob + model_order_params.is_old_const - model_order_params.is_new_const;
            else % k is new, merge_k is old
                prob = prob - model_order_params.is_old_const + model_order_params.is_new_const;
            end
        end
    end % if merge_k==0 then prob is unchanged
end