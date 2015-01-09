function neighbors = U_find_border_SP(IMG_label, IMG_SP, k, neighbors)
    [xdim, ydim] = size(IMG_label);
    found_borders = find(IMG_SP(k).borders);
    for bindex=1:length(found_borders)
        index = found_borders(bindex);
        [x, y] = get_x_and_y_from_index(index, xdim);
        if (x>1 && IMG_label(x-1, y)~=k && IMG_label(x-1, y) > 0)
            neighbors(IMG_label(x-1, y)) = true;
        end
        if (y>1 && IMG_label(x, y-1)~=k && IMG_label(x, y-1) > 0)
            neighbors(IMG_label(x, y-1)) = true;
        end
        if (x<xdim && IMG_label(x+1, y)~=k && IMG_label(x+1, y) > 0)
            neighbors(IMG_label(x+1, y)) = true;
        end
        if (y<ydim && IMG_label(x, y+1)~=k && IMG_label(x, y+1) > 0)
            neighbors(IMG_label(x, y+1)) = true;
        end
    end
end