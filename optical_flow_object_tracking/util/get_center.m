function center = get_center(mask)
    half_area = sum(mask(:))/2;
    density_column = sum(mask, 2);
    total = 0;
    row = 0;
    while total < half_area
        row = row + 1;
        total = total + density_column(row);
    end

    density_row = sum(mask, 1);
    total = 0;
    col = 0;
    while total < half_area
        col = col + 1;
        total = total + density_row(col);
    end
    center = [row col];
end