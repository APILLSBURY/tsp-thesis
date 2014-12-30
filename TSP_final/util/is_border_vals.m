% =============================================================================
% == is_border_vals.m
% == --------------------------------------------------------------------------
% == A fucntion to find the borders of the superpixels
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-4-2014
% =============================================================================

function border = is_border_vals(labels)
    [xdim, ydim] = size(labels);
    border = false(xdim, ydim);
    for x=1:xdim
        for y=1:ydim
            %make sure we don't reference outside of the matrix
            x1 = max(1, x-1);
            y1 = max(1, y-1);
            x2 = min(xdim, x+1);
            y2 = min(ydim, y+1);

            %check to see if any of the adjacent pixels are different than
            %the current one
            curLabel = labels(x, y);
            border(x, y) = (curLabel ~= labels(x1, y) || curLabel ~= labels(x, y1) || ...
                curLabel ~= labels(x2, y) || curLabel ~= labels(x, y2));
        end
    end
end