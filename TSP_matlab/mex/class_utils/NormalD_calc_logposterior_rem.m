function logprob = NormalD_calc_logposterior_rem(nd, data)
   nd = NormalD_update_stats_rem(nd, data);
   logprob = NormalD_calc_logposterior(nd, true);
end