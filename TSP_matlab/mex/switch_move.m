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
    empty_SPs = false(IMG.K, 1);
    for k=1:IMG.K
        if SP_is_empty(IMG, k) && IMG.SP_old(k)
            empty_SPs(k) = true;
        end
    end

    for k=1:IMG.K
        if SP_is_empty(IMG, k)
            % if old, check to see a new or unused old
            % if new, check to see if any unused old
            best_k = -1;
            best_energy = 0;
            if IMG.SP_old(k)
                delta = move_switch_calc_delta(IMG, IMG.SP(k), create_SP_new(IMG), true, false);
                if (delta > best_energy)
                    best_k = IMG.K+1;
                    best_energy = delta;
                end
            end

            for test_k=find(empty_SPs)'
                delta = move_switch_calc_delta(IMG, IMG.SP(k), IMG.SP(test_k), IMG.SP_old(k), IMG.SP_old(test_k));
                if (delta > best_energy)
                    best_k = test_k;
                    best_energy = delta;
                end
            end

            % switch with best label
            if best_k>0 && best_energy>0
                disp('Switchin!');
                % change the labels
                for index=find(IMG.SP(k).pixels)'
                    [x, y] = get_x_and_y_from_index(index, IMG.xdim);
                    IMG.label(x, y) = best_k;
                end
                IMG.SP_changed(k) = true;
                IMG.SP_changed(best_k) = true;

                % make room for the new one
                if best_k==IMG.K+1
                    IMG.K = IMG.K+1;
                    if IMG.K>IMG.max_SPs
                        disp('Ran out of space!');
                    end
                    if SP_is_empty(IMG, IMG.K)
                        IMG.SP(IMG.K) = new_SP(IMG.new_pos, IMG.new_app, IMG.max_UID, [0, 0], IMG.N, IMG.max_SPs);
                        IMG.max_UID = IMG.max_UID + 1;
                    end
                else
                    IMG.alive_dead_changed = true;
                end
                
                IMG = U_merge_SPs(IMG, best_k, IMG.K);
                
                % delete it if it was a new one
                if IMG.SP_old(k)
                    % add it to the search list
                    empty_SPs(k) = true;
                    IMG.alive_dead_changed = true;
                end

                %ELSE DELETE IMG.SP(k)?
                
                % remove it from the search list
                if IMG.SP_old(best_k)
                    empty_SPs(best_k) = false;
                end

                % now update the neighbors list
                % IMG = U_fix_neighbors_self(IMG, best_k);
                % update the neighbors' neighbors
                % IMG = U_fix_neighbors_neighbors(IMG, best_k, IMG.N+1);
            end
        end
    end
end


function logprob = move_switch_calc_delta(IMG, oldSP, newSP, oldSP_is_old, newSP_is_old)
    logprob = SP_log_likelihood_switch_prior(oldSP, newSP, newSP_is_old);
    logprob = logprob - oldSP.log_likelihood();
    logprob = logprob + U_calc_model_order(IMG, oldSP.N, newSP_is_old);
    logprob = logprob - U_calc_model_order(IMG, oldSP.N, oldSP_is_old);
end
