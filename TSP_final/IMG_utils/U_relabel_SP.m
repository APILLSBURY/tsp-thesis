function IMG = U_relabel_SP(IMG, final)
    last_k = 1;
    for k=1:numel(IMG.SP)
        if ~final || IMG.SP(k).N>0
            if (k~=last_k)
                % move the super pixel to the empty one
                IMG.SP(last_k) = IMG.SP(k);
                IMG.SP(k) = [];

                % relabel the pixels
                for i=1:length(IMG.SP(last_k).pixels)
                    if IMG.SP(last_k).pixels(i)
                        [x, y] = get_x_and_y_from_index(index, IMG.xdim);
                        IMG.label(x, y) = last_k;
                    end
                end
            end
            last_k = last_k + 1;
        end
    end
    IMG.K = last_k-1;
end