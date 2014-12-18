function logprob = NormalD_calc_logposterior_new2(nd, nd_other, nd_other2)
   nd = NormalD_update_stats_new2(nd, nd_other, nd_other2);
   logprob = NormalD_calc_logposterior(nd, true);
end