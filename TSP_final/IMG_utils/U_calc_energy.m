function E = U_calc_energy(IMG_N, IMG_SP, IMG_SP_old, model_order_params, IMG_new_SP)
    E = 0;
    numEmpty = 0;
    for k=1:IMG_N
        if ~(k > numel(IMG_SP) || isempty(IMG_SP(k).N) || IMG_SP(k).N == 0)
            IMG_SP(k) = SP_calculate_log_probs(IMG_SP(k));
            
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
            
            E = E + IMG_SP(k).log_likelihood + model_order;
        else
            numEmpty = numEmpty + 1;
        end
    end
    
    E = E + numEmpty*IMG_new_SP.log_likelihood;
end