% --------------------------------------------------------------------------
% -- calc_logposterior_MM
% --   calculates log p(x | params)
% --------------------------------------------------------------------------
function logprob = NormalD_calc_logposterior_MM(nd, data)
    mean = NormalD_calc_mean(nd);
    logprob = -nd.D*0.5*1.837877066409 - ((data - mean).^2)/2;
end