function model_order = U_calc_model_order(IMG_log_beta, IMG_log_area_var, IMG_area, IMG_area_var, IMG_log_alpha, size, is_old)
    if size==0
        model_order = 0;
    elseif is_old
        model_order = IMG_log_beta - 0.5  * (1.837877066409 + IMG_log_area_var + (IMG_area - size)^2/IMG_area_var);
    else
        model_order = IMG_log_alpha - 0.5 * (1.837877066409 + IMG_log_area_var + (IMG_area - size)^2/IMG_area_var);
    end
    % model_order_params.is_old/new_const + (model_order_params.area - size/2)*size/model_order_params.area_var

    %calculate model order
    size = IMG_SP(k).N;
    if size==0
        model_order = 0;
    else
        size_var = (model_order_params.area - size/2)*size/model_order_params.area_var;
        if IMG_SP_old(k)
            model_order = model_order_params.is_old_const + size_var;
        else
            model_order = model_order_params.is_new_const + size_var;
        end
    end
end

