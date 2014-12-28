% =============================================================================
% == merge_move.m
% == --------------------------------------------------------------------------
% == Finds the optimal merge between all pairs of super pixels.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-4-2014
% =============================================================================

function IMG = merge_move(IMG, its)
    addpath('mex/class_utils/');

    for i=1:its
        % choose a random order of super pixels
        Nsp = numel(IMG.SP);
        perm = randperm(Nsp);

        neighbors = false(Nsp, 1);

        for k=perm;

            % find a nonempty super pixel
            if ~SP_is_empty(IMG, k)
                % find all bordering super pixels
                neighbors = U_find_border_SP(IMG, k, neighbors);

                max_E = -inf;
                max_k = -1;

                % loop through all neighbors
                for merge_k=find(neighbors)'
                    new_E = move_merge_calc_delta(IMG, k, merge_k);
                    if new_E > max_E || max_k==-1
                        max_E = new_E;
                        max_k = merge_k;
                    end
                end

                % fix neighbors for the next super pixel
                neighbors = false(size(neighbors));

                % merge if it increases energy
                if (max_E>0)
                    disp('max e is greater than 0');
                    % change the labels
                    for index=1:length(IMG.SP(max_k).pixels)
                        if IMG.SP(max_k).pixels(index)
                            [x, y] = get_x_and_y_from_index(index, IMG.xdim);
                            IMG.label(x, y) = k;
                        end
                    end
                    IMG.SP_changed(k) = true;
                    IMG.SP_changed(max_k) = true;

                    IMG.SP(k) = SP_merge_with(IMG.SP(k), IMG.SP(max_k), IMG.label);
                    IMG.SP(max_k) = SP_empty(IMG.SP(max_k));
                    if (~IMG.SP_old(max_k))
                        IMG.SP(max_k) = SP_empty(IMG.SP(max_k)); %VERY QUESTIONABLE
                    else
                        IMG.alive_dead_changed = true;
                    end
                end
            end
        end
    end
end

