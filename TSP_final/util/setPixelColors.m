function [output] = setPixelColors(input, indices, color)
    addpath('mex/class_utils/');
    [xdim, ydim] = size(input);
    output = input;
    for index=indices'
        [x, y] = get_x_and_y_from_index(index, xdim);
        output(x, y, 1) = color(1);
        output(x, y, 2) = color(2);
        output(x, y, 3) = color(3);
    end
end