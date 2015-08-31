% implements a depth first search to find all connected areas in a boolean
% matrix
function object_masks = find_connected_areas(boolean_mask)
    [xdim, ydim] = size(boolean_mask);
    object_masks = logical.empty(xdim, ydim, 0);
    x=1;
    level = 0;
    while x<=xdim
        temp_y = find(boolean_mask(x,:), 1);
        if isempty(temp_y)
            x = x+1;
            continue;
        end
        temp_x = x;
        bfs_write = 1;
        bfs_read = 1;
        bfs = zeros(xdim*ydim, 2);
        bfs(bfs_write, :) = [temp_x, temp_y];
        level = level + 1;
        object_masks(:, :, level) = false(xdim, ydim); % create a new level
        object_masks(temp_x, temp_y, level) = true;
        boolean_mask(temp_x, temp_y) = false;
        bfs_write = bfs_write+1;
        while bfs_read < bfs_write
            temp_x = bfs(bfs_read, 1);
            temp_y = bfs(bfs_read, 2);
            if (temp_x>1 && boolean_mask(temp_x-1, temp_y))
                bfs(bfs_write, :) = [temp_x-1, temp_y];
                object_masks(temp_x-1, temp_y, level) = true;
                boolean_mask(temp_x-1, temp_y) = false;
                bfs_write = bfs_write+1;
            end
            if (temp_y>1 && boolean_mask(temp_x, temp_y-1))
                bfs(bfs_write, :) = [temp_x, temp_y-1];
                object_masks(temp_x, temp_y-1, level) = true;
                boolean_mask(temp_x, temp_y-1) = false;
                bfs_write = bfs_write+1;
            end
            if (temp_x<xdim && boolean_mask(temp_x+1, temp_y))
                bfs(bfs_write, :) = [temp_x+1, temp_y];
                object_masks(temp_x+1, temp_y, level) = true;
                boolean_mask(temp_x+1, temp_y) = false;
                bfs_write = bfs_write+1;
            end
            if (temp_y<ydim && boolean_mask(temp_x, temp_y+1))
                bfs(bfs_write, :) = [temp_x, temp_y+1];
                object_masks(temp_x, temp_y+1, level) = true;
                boolean_mask(temp_x, temp_y+1) = false;
                bfs_write = bfs_write+1;
            end
            bfs_read = bfs_read+1;
        end
    end
end