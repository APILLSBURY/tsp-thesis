% finds the maximum overlap between a moving mask and a stationary mask
% given the direction that the moving mask is going
% USAGE: theta_median should be the theta for whichever is the later mask
% set reverse_theta to true when stationary mask is the later mask
function [max_offset_mask, max_row_offset, max_col_offset] = get_overlap(theta_median, moving_mask, stationary_mask, reverse_theta)
    if reverse_theta
        theta_median = theta_median + pi;
        if theta_median>pi
            theta_median = theta_median - 2*pi;
        end
    end

    [xdim, ydim] = size(moving_mask);
    max_overlap = -1;
    max_offset_mask = false(xdim, ydim);
    max_row_offset = 0;
    max_col_offset = 0;
    rho = 0;
    found_overlap = false;
    while rho==0 || sum(offset_mask(:))>0 && ~(overlap==0 && found_overlap)
        [col_offset, row_offset] = pol2cart(theta_median, rho);
        row_offset = round(row_offset);
        col_offset = round(col_offset);
        offset_mask = get_offset_mask(moving_mask, row_offset, col_offset);
        combined_mask = offset_mask & stationary_mask;
        overlap = sum(combined_mask(:));
        if overlap>0
            found_overlap = true;
        end
        if overlap>max_overlap
            max_overlap = overlap;
            max_offset_mask = offset_mask;
            max_row_offset = row_offset;
            max_col_offset = col_offset;
        end
        rho = rho + 5;
    end
end