function IMG = IMG_U_initialize(IMG)
    IMG.K = max(max(IMG.label));
    IMG.new_pos = new_NormalD(2, IMG.hyper.p_theta, IMG.hyper.p_Delta, true);
    IMG.new_app = new_NormalD(3, IMG.hyper.a_theta, IMG.hyper.a_Delta, true);
    IMG.SP_old = false(1,IMG.N);
    IMG.log_area_var = log(IMG.area_var);
    for i=1:IMG.K
        IMG.SP(i) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, [0, 0], IMG.N);
        IMG.max_UID = IMG.max_UID+1;
    end

    % populate the linked lists and pointers
    for x=1:IMG.xdim
        for y=1:IMG.ydim
            curLabel = IMG.label(x, y);
            if (curLabel>0)
                index = get_index_from_x_and_y(x, y, IMG.xdim);
                IMG.SP(curLabel) = SP_add_pixel_init(IMG.SP(curLabel), IMG.data, index, U_check_border_pix(IMG, index), IMG.boundary_mask(x, y));
                IMG.SP(curLabel) = SP_update_neighbors_add_self(IMG.SP(curLabel), IMG.label, index);
            end
        end
    end
    
    for k=1:IMG.K
        IMG.SP(k) = SP_calculate_log_probs(IMG.SP(k));
    end
end