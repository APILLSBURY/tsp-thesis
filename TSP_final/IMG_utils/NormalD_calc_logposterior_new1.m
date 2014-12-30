% --------------------------------------------------------------------------
% -- calc_logposterior_new
% --   calculate log p(params | x, hyperparams), for the maximal params
% -- with a new merge. Because of conjugacy, this is in the same class
% -- as log p(params | hyperparams)
% --
% --   parameters:
% --     - other : another NormalD to merge with
% --------------------------------------------------------------------------
function logprob = NormalD_calc_logposterior_new1(nd, nd_other)
   nd = NormalD_update_stats_new1(nd, nd_other);
   logprob = NormalD_calc_logposterior(nd, true);
end