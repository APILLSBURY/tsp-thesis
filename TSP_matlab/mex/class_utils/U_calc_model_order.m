function model_order = U_calc_model_order(IMG, size, is_old)
    if size==0
        model_order = 0;
    elseif is_old
        model_order = IMG.log_beta - 0.5*1.837877066409 - 0.5*IMG.log_area_var - 0.5*(IMG.area-size)^2/IMG.area_var;
    else
        model_order = IMG.log_alpha - 0.5*1.837877066409 - 0.5*IMG.log_area_var - 0.5*(IMG.area-size)^2/IMG.area_var;
    end
end