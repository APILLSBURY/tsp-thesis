% --------------------------------------------------------------------------
% -- calc_logposterior_new
% --   calculate log p(params | x, hyperparams), for the maximal params
% -- with a new data point. Because of conjugacy, this is in the same class
% -- as log p(params | hyperparams)
% --
% --   parameters:
% --     - data : the new point to possibly add
% --------------------------------------------------------------------------
function logprob = NormalD_calc_logposterior_new(nd, data)
   nd = NormalD_update_stats_new(nd, data);
   logprob = NormalD_calc_logposterior(nd, true);
end