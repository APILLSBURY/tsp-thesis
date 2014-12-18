% =============================================================================
% == switch_move.cpp
% == --------------------------------------------------------------------------
% == A MEX interface to perform switch moves on TSPs.
% == See m files for calling convention.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-12-2014
% =============================================================================

function IMG = switch_move(IMG)
    addpath('mex/class_utils/');
    empty_SPs = false(length(IMG.SP), 1);
    for k=1:length(IMG.SP)
        if (isempty(IMG.SP(k).N) && IMG.SP_old(k))
            empty_SPs(k) = true;
        end
    end

    for k=1:length(IMG.SP)
        if (~isempty(IMG.SP(k).N))
            % if old, check to see a new or unused old
            % if new, check to see if any unused old
            best_k = -1;
            best_energy = 0;
            if (IMG.SP_old(k))
                delta = move_switch_calc_delta(IMG.SP(k), create_SP_new(IMG));
                if (delta > best_energy)
                    best_k = IMG.K;
                    best_energy = delta;
                end
            end

            for test_k=1:length(empty_SPs)
                if empty_SPs(test_k)
                    delta = move_switch_calc_delta(IMG.SP(k), IMG.SP(test_k));
                end
                if (delta > best_energy)
                    best_k = test_k;
                    best_energy = delta;
                end
            end

            % switch with best label
            if (best_k>0 && best_energy>0)
                % change the labels
                for index = 1:length(IMG.SP(k).pixels)
                    if IMG.SP(k).pixels(index)
                        [x, y] = get_x_and_y_from_index(index, IMG.xdim);
                        IMG.label(x, y) = best_k;
                    end
                end
                IMG.SP_changed(k) = true;
                IMG.SP_changed(best_k) = true;

                % make room for the new one
                if (best_k==IMG.K)
                    if (IMG.K>=IMG.N)
                        disp('Ran out of space!');
                    end
                    if (isempty(IMG.SP(k).N))
                        IMG.SP(k) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, false, IMG.SP(k).prev_v, IMG.N);
                        IMG.max_UID = IMG.max_UID + 1;
                    end
                    IMG.K = IMG.K + 1;
                else
                    IMG.alive_dead_changed = true;
                end
                IMG.SP(best_k) = SP_merge_with(IMG.SP(best_k), IMG.SP(k), IMG.label);
                IMG.SP(k) = IMG.SP_empty(IMG.SP(k));
                
                % delete it if it was a new one
                if (~IMG.SP_old(k))
                    IMG.SP(k) = [];
                else
                    % add it to the search list
                    empty_SPs(k) = true;
                    IMG.alive_dead_changed = true;
                end

                % remove it from the search list
                if (IMG.SP_old(best_k))
                    empty_SPs(best_k) = false;
                end

                % now update the neighbors list
                IMG = U_fix_neighbors_self(IMG, best_k);
                % update the neighbors' neighbors
                IMG = U_fix_neighbors_neighbors(IMG, best_k, IMG.N+1);
            end
        end
    end
end


function logprob = move_switch_calc_delta(IMG, oldSP, newSP)
    logprob = SP_log_likelihood_switch_prior(oldSP, newSP);
    logprob = logprob - oldSP.log_likelihood();
    logprob = logprob + U_calc_model_order(IMG, oldSP.N, oldSP.is_old);
    logprob = logprob - U_calc_model_order(IMG, oldSP.N, newSP.is_old);
end
