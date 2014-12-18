% --------------------------------------------------------------------------
% -- calc_logposterior_MM
% --   calculates log p(x | params)
% --------------------------------------------------------------------------
function logprob = NormalD_calc_logposterior_MM(nd, data)
    logprob = -nd.D*0.5*1.837877066409 - ((data - nd.mean).^2)/2;
end