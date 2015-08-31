% finds the mask that is offset by a given number of rows and columns
function offset_mask = get_offset_mask(moving_mask, row_offset, col_offset)
    offset_mask = false(size(moving_mask));
    offset_mask(max(1+row_offset, 1):min(end, end+row_offset), max(1+col_offset, 1):min(end, end+col_offset)) = moving_mask(max(1, 1-row_offset):min(end-row_offset, end), max(1, 1-col_offset):min(end-col_offset, end));
end