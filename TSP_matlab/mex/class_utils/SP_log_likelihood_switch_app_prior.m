% --------------------------------------------------------------------------
% -- log_likelihood_switch_app_prior
% --   Calculates the likelihood for this SP if the appearance prior is
% -- switched to the prior of the supplied SP
% --
% --   parameters:
% --     - new_prior: a pointer to the SP to switch priors with
% --------------------------------------------------------------------------
function logprob = SP_log_likelihood_switch_app_prior(SP, SP_new_prior)
   total = SP.app.total;
   total2 = SP.app.total2;
   logprob = -3*(SP.N+0.5)*0.5*1.837877066409 - SP_new_prior.app.sumlogDelta_div2;
   if ~SP_new_prior.is_old
       logprob = logprob + sum(total.*total./SP.N - total2)/2;
   else
       Delta = SP_new_prior.app.Delta;
       theta = SP_new_prior.app.theta;
       logprob = logprob + sum((Delta.*total.*total + 2.*total.*theta - SP.N.*theta.*theta) ./ (2.*SP.N.*Delta+2) - total2./2);
   end
end