function [IMG_SP, IMG_label, IMG_K] = U_relabel_SP(IMG_SP, IMG_label, final)
    [xdim, ~] = size(IMG_label);
    last_k = 1;
    for k=1:numel(IMG_SP)
        if ~final || IMG_SP(k).N>0
            if (k~=last_k)
                % move the super pixel to the empty one
                IMG_SP(last_k) = IMG_SP(k);
                IMG_SP(k) = [];

                % relabel the pixels
                found_pix = find(IMG_SP(last_k).pixels);
                for index=1:length(found_pix)
                    [x, y] = get_x_and_y_from_index(found_pix(index), xdim);
                    IMG_label(x, y) = last_k;
                end
            end
            last_k = last_k + 1;
        end
    end
    IMG_K = last_k-1;
end