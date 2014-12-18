% --------------------------------------------------------------------------
% -- log_likelihood_switch_prior
% --   Calculates the likelihood for this SP if the prior is switched to
% -- the prior of the supplied SP
% --
% --   parameters:
% --     - new_prior: a pointer to the SP to switch priors with
% --------------------------------------------------------------------------
function prob = SP_log_likelihood_switch_prior(SP, SP_new_prior)
   prob = SP_log_likelihood_switch_pos_prior(SP, SP_new_prior) + SP_log_likelihood_switch_app_prior(SP, SP_new_prior);
end