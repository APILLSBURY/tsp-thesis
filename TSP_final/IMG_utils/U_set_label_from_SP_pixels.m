function IMG_label = U_set_label_from_SP_pixels(IMG_label, SP, curlabel)
    found_pixels = find(SP.pixels);
    [xdim, ~] = size(IMG_label);
    for index=1:length(found_pixels)
        pix = found_pixels(index);
        [x, y] = get_x_and_y_from_index(pix, xdim);
        IMG_label(x, y) = curlabel;
    end
end