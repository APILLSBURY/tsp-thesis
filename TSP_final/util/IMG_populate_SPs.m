function IMG = IMG_populate_SPs(IMG)
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