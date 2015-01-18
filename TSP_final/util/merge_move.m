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

function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_alive_dead_changed, IMG_SP_old]  = merge_move(IMG_label, IMG_SP, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, model_order_params, IMG_K, its)
    xdim = size(IMG_label, 1);
    for i=1:its
        % choose a random order of super pixels
        Nsp = numel(IMG_SP);
        perm = randperm(Nsp);

        neighbors = false(Nsp, 1);

        for kindex=1:length(perm)
            k=perm(kindex);

            % find a nonempty super pixel
            if ~(k > numel(IMG_SP) || isempty(IMG_SP(k).N) || IMG_SP(k).N == 0)
                % find all bordering super pixels
                neighbors = U_find_border_SP(IMG_label, IMG_SP, k, neighbors);

                max_E = -inf;
                max_k = -1;

                % loop through all neighbors
                found_neighbors = find(neighbors);
                for merge_k_index=1:length(found_neighbors);
                    merge_k = found_neighbors(merge_k_index);
                    new_E = U_move_merge_calc_delta(IMG_SP(k), IMG_SP(merge_k), IMG_SP_old(k), IMG_SP_old(merge_k), model_order_params);
                    if new_E > max_E || max_k==-1
                        max_E = new_E;
                        max_k = merge_k;
                    end
                end

                % fix neighbors for the next super pixel
                neighbors = false(size(neighbors));

                fprintf('max_E = %f\n', max_E);
                % merge if it increases energy
                if (max_E>0)
                    disp('max e is greater than 0');
                    % change the labels
                    found_pixels = find(IMG_SP(max_k).pixels);
                    for index=1:found_pixels
                        [x, y] = get_x_and_y_from_index(found_pixels(index), xdim);
                        IMG_label(x, y) = k;
                    end
                    IMG_SP_changed(k) = true;
                    IMG_SP_changed(max_k) = true;

                    IMG_SP = U_merge_SPs(IMG_SP, IMG_label, k, max_k);
                    if IMG_SP_old(max_k)
                        IMG_alive_dead_changed = true;
                    elseif max_k<numel(IMG_SP)
                        IMG_label = U_set_label_from_SP_pixels(IMG_label, IMG_SP(numel(IMG_SP)), max_k);
                        IMG_SP = U_merge_SPs(IMG_SP, IMG_label, max_k, numel(IMG_SP));
                        IMG_SP(numel(IMG_SP)) = [];
                        IMG_K = IMG_K-1;
                    end
                end
            end
        end
    end
end

