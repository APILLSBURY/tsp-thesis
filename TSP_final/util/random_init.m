function label = random_init(xdim, ydim, w, K)
    disp('starting random_init');
    rxdim = xdim - 2*w;
    rydim = ydim - 2*w;

    N = rxdim*rydim;

    centers = randsample(N, K);
    %centers = [32630; 20589; 100; 30251; 10943; 18707; 14529; 54253; ...
    %            22853; 13432; 16948; 261; 5691; 9156; 38022; 60294];
    
    [centers_x, centers_y] = ind2sub([rxdim,rydim], centers);
    centers = [centers_x centers_y] + w;
    label = zeros(xdim, ydim);

    %first column is x, second is y, third is label
    bfs_queue = zeros(xdim*ydim*4, 3);
    bfs_counter = 1;

    %set the label centers
    for k=1:K
        label(centers(k,1), centers(k,2)) = k;
        bfs_queue(bfs_counter, :) = [centers(k,1), centers(k,2), k];
        bfs_counter = bfs_counter+1;
    end

    bfs_index = 1;
    while bfs_index < bfs_counter
        x = bfs_queue(bfs_index, 1);
        y = bfs_queue(bfs_index, 2);
        k = bfs_queue(bfs_index, 3);

        if x>1 && label(x-1, y)==0
            label(x-1, y) = k;
            bfs_queue(bfs_counter, :) = [x-1, y, k];
            bfs_counter = bfs_counter+1;
        end
        if y>1 && label(x, y-1)==0
            label(x, y-1) = k;
            bfs_queue(bfs_counter, :) = [x, y-1, k];
            bfs_counter = bfs_counter+1;
        end
        if x<xdim && label(x+1, y)==0
            label(x+1, y) = k;
            bfs_queue(bfs_counter, :) = [x+1, y, k];
            bfs_counter = bfs_counter+1;
        end
        if y<ydim && label(x, y+1)==0
            label(x, y+1) = k;
            bfs_queue(bfs_counter, :) = [x, y+1, k];
            bfs_counter = bfs_counter+1;
        end
        bfs_index = bfs_index+1;
    end
end
