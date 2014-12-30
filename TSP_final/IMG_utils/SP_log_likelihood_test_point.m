% --------------------------------------------------------------------------
% -- log_likelihood_test_pointXXXX
% --   finds the log likelihood for adding or removing a data point. All
% -- functions take one parameter, a 5-d double vector of data. Some
% -- functions indicated by * take an additional boolean that indicates if
% -- the position should be checked, or the position and the appearance.
% -- Defaults to true, which checks both.
% -- 
% -- XXXX can be any of the following
% --
% -- <empty>* : calculate the likelihood for adding
% -- rem* : calculate the likelihood for removing
% -- app : calculates the appearance likelihood for adding
% -- pos : calculates the position likelihood for adding
% -- app_rem : calculates the appearance likelihood for removing
% -- pos_rem : calculates the position likelihood for removing
% -- MM : appearance and position likelihood for adding w/out finding the
% --   corresponding optimal parameters, similar to an iterative scheme
% -- MM_pos : same as above, except only for position
% -- 
% --   parameters
% --     - data : a [1 5] vector of a test point to add
% --------------------------------------------------------------------------
function prob = SP_log_likelihood_test_point(SP, data, checkApp)
   % calculate the new log probability
   if (checkApp)
      prob = NormalD_calc_logposterior_new(SP.pos, data(1:2)) + NormalD_calc_logposterior_new(SP.app, data(3:5));
   else
      prob = NormalD_calc_logposterior_new(SP.pos, data(1:2)) + NormalD_calc_logposterior(SP.app);
   end
end