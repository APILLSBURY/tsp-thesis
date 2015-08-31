% expand a mask in every direction by the expansion variable
function expanded_mask = expand_mask(mask, expansion)
    [xdim, ydim] = size(mask);
    expanded_mask = mask;
    %expand it horizontally
    for x=1:xdim
        found_row = find(mask(x,:));
        prev = 0;
        first = 0;
        for y_index=1:length(found_row)
            y = found_row(y_index);
            if prev == 0
                first = y;
            elseif y-prev~=1 %we found a gap
                expanded_mask(x, max(1, first-expansion):min(ydim, prev+expansion)) = true;
                first = y;
            end
            prev = y;
        end
        %finish off the last one
        if first~=0
            expanded_mask(x, max(1, first-expansion):min(ydim, prev+expansion)) = true;
        end
    end

    %expand it vertically
    for y=1:ydim
        found_column = find(mask(:,y));
        prev = 0;
        first = 0;
        for x_index=1:length(found_column)
            x = found_column(x_index);
            if prev == 0
                first = x;
            elseif x-prev~=1 %we found a gap
                expanded_mask(max(1, first-expansion):min(xdim, prev+expansion), y) = true;
                first = x;
            end
            prev = x;
        end
        %finish off the last one
        if first~=0
            expanded_mask(max(1, first-expansion):min(xdim, prev+expansion), y) = true;
        end
    end
end