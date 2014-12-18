% --------------------------------------------------------------------------
% -- add_data
% --   functions to add an observation to the NormalD. Updates the sufficient
% -- statistics stored in the data structure, but not the posterior
% -- hyperparameters
% --
% --   parameters:
% --     - data : the new observed data point of size [1 D]
% --------------------------------------------------------------------------
function nd = NormalD_add_data(nd, data)
    nd.N = nd.N + 1;
    nd.total = nd.total + data;
    nd.total2 = nd.total2 + data.^2;
end