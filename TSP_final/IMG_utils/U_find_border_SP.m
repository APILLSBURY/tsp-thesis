function neighbors = U_find_border_SP(IMG, k, neighbors)
    for index=find(IMG.SP(k).borders)'
        [x, y] = get_x_and_y_from_index(index, IMG.xdim);
        if (x>1 && IMG.label(x-1, y)~=k && IMG.label(x-1, y) > 0)
            neighbors(IMG.label(x-1, y)) = true;
        end
        if (y>1 && IMG.label(x, y-1)~=k && IMG.label(x, y-1) > 0)
            neighbors(IMG.label(x, y-1)) = true;
        end
        if (x<IMG.xdim && IMG.label(x+1, y)~=k && IMG.label(x+1, y) > 0)
            neighbors(IMG.label(x+1, y)) = true;
        end
        if (y<IMG.ydim && IMG.label(x, y+1)~=k && IMG.label(x, y+1) > 0)
            neighbors(IMG.label(x, y+1)) = true;
        end
    end
end