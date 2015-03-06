function [x, y] = get_x_and_y_from_index(index, xdim)
    x = mod(index,xdim);
    if x==0
        x=xdim;
    end
    y = ceil(index/xdim);
end