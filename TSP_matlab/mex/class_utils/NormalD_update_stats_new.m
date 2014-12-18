% --------------------------------------------------------------------------
% -- update_stats_new
% --   updates the internal posterior hyperparameters based on the stored
% -- sufficient statistics and the new data points
% --
% --   parameters:
% --     - data: the new point to consider
% --------------------------------------------------------------------------
function nd = NormalD_update_stats_new(nd, data)
   nd.temp_N = nd.N + 1;
   nd.temp_total = nd.total + data;
   nd.temp_total2 = nd.total2 + data.^2;
end