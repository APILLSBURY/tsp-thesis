function IMG = U_set_label_from_SP_pixels(IMG, SP, label)
    for pix=find(SP.pixels)'
        [x, y] = get_x_and_y_from_index(pix, IMG.xdim);
        IMG.label(x, y) = label;
    end
end