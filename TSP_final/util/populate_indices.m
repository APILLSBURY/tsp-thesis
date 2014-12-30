% =============================================================================
% == calc_sdfIMPORT.cpp
% == --------------------------------------------------------------------------
% == The MEX interface file to calculate a signed distance function
% == --------------------------------------------------------------------------
% == Copyright 2011. MIT. All Rights Reserved.
% == Written by Jason Chang 06-13-2011
% =============================================================================

function indices = populate_indices(K, label)
    addpath('util/');
    
    [xdim, ydim] = size(label);
    indices = repmat(struct('all', zeros(1, xdim*ydim)), K, 1);
    counter = zeros(K,1);
    for x=1:xdim
        for y=1:ydim
            index = get_index_from_x_and_y(x, y, xdim);
            k = label(x, y);
            if (k>0)
                counter(k)=counter(k)+1;
                indices(k).all(counter(k)) = index;
            end
        end
    end
    for k=1:K
        indices(k).all = indices(k).all(1:counter(k));
    end
end
