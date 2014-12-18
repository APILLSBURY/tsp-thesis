% --------------------------------------------------------------------------
% -- log_likelihood_test_merge
% --   Calculate the log likelihood for merging with another SP. Can be 
% -- called with two SPs as arguments, in which case, it checks the merge
% -- of all three SPs.  Can also use _pos to check only the position terms.
% --
% --   parameters
% --     - other : another SP that we are testing for a merge
% --------------------------------------------------------------------------
function prob = SP_log_likelihood_test_merge1(SP, SP_other)
    % calculate the new log probability
    prob = NormalD_calc_logposterior_new1(SP.pos, SP_other.pos) + NormalD_calc_logposterior_new1(SP.app, SP_other.app);
end
