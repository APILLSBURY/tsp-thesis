function prob = SP_log_likelihood_test_point_rem(SP, data, checkApp)
   % calculate the new log probability
   if (checkApp)
       prob = NormalD_calc_logposterior_rem(SP.pos, data(1:2)) + NormalD_calc_logposterior_rem(SP.app, data(3:5));
   else
      prob = NormalD_calc_logposterior_new(SP.pos, data(1:2)) + NormalD_calc_logposterior(SP.app);
   end
end