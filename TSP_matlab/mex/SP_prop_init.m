% =============================================================================
% == SP_prop_init.cpp
% == --------------------------------------------------------------------------
% == A MEX interface to propagate labels according to some flow.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-12-2014
% =============================================================================

function label = SP_prop_init(K, label, meanx, meany, vx, vy, mask)

    [xdim, ydim] = size(label);
    N = xdim*ydim;

    newlabel = zeros(xdim, ydim);
    possible_labels = false(xdim, ydim, K);
    for x=1:xdim
        for y=1:ydim
            %initialize newLabel
            if mask(x, y)
                newlabel = K;
            end

            %find possible labels
            if label(x, y) > 0
                k = label(x, y);
                newx = x + vx(k) + 0.5;
                newy = y + vy(k) + 0.5;

                if newx>0 && newx<xdim && newy>0 && newy<ydim
                    possible_labels(x, y, k) = true;
                end
            end
        end
    end

    %propagate labels
    newRegion = false;
    for x=1:xdim
        for y=1:ydim
            possibilities = find(possible_labels(x, y, :));
            if length(possibilities)==1
                newlabel(x, y) = possibilities(1);
            elseif length(possibilities)>1
                closest_k = -1;
                closest_distance = N*N;
                for k=possibilities
                    distance = (meanx(k) + vx(k) - x)^2 + (meany(k) + vy(k) - y)^2;
                    if (distance<closest_distance || closest_k<0)
                        closest_distance = distance;
                        closest_k = k;
                    end
                end
                newlabel(x, y) = closest_k;
            elseif mask(x, y)
                newRegion = true;
            end
        end
    end
    if (newRegion)
        K = K+1;
    end

    % labels propagated, now enforce connectivity
    done = false(xdim, ydim);
    SP_to_group_label = zeros(K, N);
    group_label = zeros(xdim, ydim);
    group_label_index = 1;
    for x=1:xdim
        for y=1:ydim
            if ~done(x, y) && newlabel(x, y)>0
                done(x, y) = true;
                curLabel = newlabel(x, y);

                explore = zeros(N, 2);
                explore_write_index = 1;
                explore(explore_write_index,:) = [x, y];
                explore_write_index = explore_write_index+1;

                explore_read_index = 1;
                while explore_read_index < explore_write_index

                    xi = explore(explore_read_index, 1);
                    yi = explore(explore_read_index, 2);

                    group_label(xi, yi) = group_label_index;

                    if (xi>1 && ~done(x-1, y) && newlabel(x-1, y)==curLabel)
                        explore(explore_write_index, :) = [x-1, y];
                        explore_write_index = explore_write_index+1;
                        done(x-1, y) = true;
                    end
                    if (yi>1 && ~done(x, y-1) && newlabel(x, y-1)==curLabel)
                        explore(explore_write_index, :) = [x, y-1];
                        explore_write_index = explore_write_index+1;
                        done(x, y-1) = true;
                    end
                    if (xi<xdim && ~done(x+1, y) && newlabel(x+1, y)==curLabel)
                        explore(explore_write_index, :) = [x+1, y];
                        explore_write_index = explore_write_index+1;
                        done(x+1, y) = true;
                    end
                    if (yi<ydim && ~done(x, y+1) && newlabel(x, y+1)==curLabel)
                        explore(explore_write_index, :) = [x, y+1];
                        explore_write_index = explore_write_index+1;
                        done(x, y+1) = true;
                    end
                    explore_read_index = explore_read_index + 1;
                end

                SP_to_group_label(curLabel, group_label_index) = explore_read_index-1;

            elseif (~done(x, y) && newlabel(x, y)<=0)
                done(x, y) = true;
            end
        end
    end

    curK = K;
    for k=1:curK
        %indices_k holds the indices of all the groups that have label k
        indices_k = find(SP_to_group_label(k, :));
        biggest = max(SP_to_group_label(k, :));
        for group=indices_k
            group_size = SP_to_group_label(k, group);
            if group_size ~= biggest && group_size > 1
                group_indices = find(group_label==group);
                if group_size<=20
                    done(group_indices) = false;
                else
                    done(group_indices) = true;
                    newlabel(group_indices) = K+1;
                    K = K+1;
                end
            end
        end
    end
    any_notdone = true;
    count = 1;
    while any_notdone && count<=10
        any_notdone = false;
        % we can either just set to a neighboring K, or create new ones
        % this sets to a neighboring k
        for x=1:xdim
            for y=1:ydim
                if ~done(x, y)
                    done(x, y) = true;
                    if (x<xdim && done(x+1,y) && newlabel(x+1,y)>0)
                        newlabel(x, y) = newlabel(x+1, y);
                    elseif (y<ydim && done(x,y+1) && newlabel(x,y+1)>0)
                        newlabel(x, y) = newlabel(x, y+1);
                    elseif (x>1 && done(x-1,y) && newlabel(x-1,y)>0)
                        newlabel(x, y) = newlabel(x-1, y);
                    elseif (y>1 && done(x,y-1) && newlabel(x,y-1)>0)
                        newlabel(x, y) = newlabel(x, y-1);
                    else
                        done(x, y) = false;
                    end
                end
            end
        end

        for x=xdim:-1:1
            for y=ydim:-1:1
                if ~done(x, y)
                    done(x, y) = true;
                    if (x<xdim && done(x+1,y) && newlabel(x+1,y)>0)
                        newlabel(x, y) = newlabel(x+1, y);
                    elseif (y<ydim && done(x,y+1) && newlabel(x,y+1)>0)
                        newlabel(x, y) = newlabel(x, y+1);
                    elseif (x>1 && done(x-1,y) && newlabel(x-1,y)>0)
                        newlabel(x, y) = newlabel(x-1, y);
                    elseif (y>1 && done(x,y-1) && newlabel(x,y-1)>0)
                        newlabel(x, y) = newlabel(x, y-1);
                    else
                        done(x, y) = false;
                        any_notdone = true;
                    end
                end
            end
        end
        count=count+1;
    end

    for x=1:xdim
        for y=1:ydim
            if ~done(x, y)
                if mask(x, y)
                    disp('SP_prop_init exceeded tries');
                else
                    newlabel(x, y) = 0;
                end
            end
        end
    end
end