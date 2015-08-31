% find how likely it is that two masks are the same object
function rating = sameness(mask, prev_mask, theta)
    % find the best overlap as a percentage of area of the larger mask
    [prev_mask_offset, ~, ~] = get_overlap(theta, prev_mask, mask, false);
    rating = sum( sum( mask & prev_mask_offset ) ) / max(sum( mask(:) ), sum( prev_mask(:) ));
end