% =============================================================================
% == IMG.m
% == --------------------------------------------------------------------------
% == An image class to be used with temporal superpixels.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-10-2014
% =============================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Utility Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Moving Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------------------------------------------------------
% -- link_cost
% --   calculates the energy of linking IMG.SP(ko) (old) to IMG.SP(kn) (new)
% --
% --   parameters:
% --     - ko : the old k
% --     - kn : the new k
% --------------------------------------------------------------------------
function logprob = link_cost(IMG, ko, kn)
    if (ko<IMG.prev_K && kn<IMG.K) % matched SPs
        logprob = SP_log_likelihood_switch_app_prior(IMG.SP(kn), IMG.SP(ko));
        logprob = logprob + U_calc_model_order(IMG, IMG.SP(kn).N, true);
    elseif (ko>=IMG.prev_K) % new SP
        logprob = SP_log_likelihood_switch_app_prior(IMG.SP(kn), create_SP_new(IMG));
        logprob = logprob + U_calc_model_order(IMG, IMG.SP(kn).N, false);
    elseif (kn>=IMG.K) % dead SP
        logprob = SP_log_likelihood_switch_app_prior(create_SP_new(IMG), IMG.SP(ko));
    else
        disp('Cant link a dead SP to a new SP');
    end
end



% --------------------------------------------------------------------------
% -- move_local_calc_neigbor
% --   calculates the probability of assigning the pixel at index to the
% -- cluster of nindex (neighbor index). Updates the max_prob and max_k if
% -- the new_app value is greater than the old one
% --
% --   parameters:
% --     - index : the new point to add
% --     - nindex : the neighbor to add this point to
% --     - max_prob : the maximum probability of all neighbors
% --     - max_k : the super pixel index of the maximum probability
% --------------------------------------------------------------------------
function [max_prob, max_k] = move_local_calc_delta_MM(IMG, index, new_k, max_prob, max_k)
    prob = SP_log_likelihood_test_point_MM(IMG.SP(new_k), IMG.data(index, :));

    if (prob>max_prob+1e-10 || max_k<0)
        max_prob = prob;
        max_k = new_k;
    end
end



