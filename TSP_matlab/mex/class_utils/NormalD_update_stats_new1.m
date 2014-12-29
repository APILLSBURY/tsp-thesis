% --------------------------------------------------------------------------
% -- update_stats_new
% --   updates the internal posterior hyperparameters based on the stored
% -- sufficient statistics and the new data points
% --
% --   parameters:
% --     - other : another NormalD to merge with and test the likelihood
% --------------------------------------------------------------------------
function nd = NormalD_update_stats_new1(nd, nd_other)
   nd.temp_N = nd.N + nd_other.N;
   nd.temp_total = nd.total + nd_other.total;
   nd.temp_total2 = nd.total2 + nd_other.total2;
end