% --------------------------------------------------------------------------
% -- NormalD
% --   constructor; intializes to all the values given
% --------------------------------------------------------------------------
function nd = new_NormalD(newD, newTheta, newDelta, newUniform)

    nd = NormalD_initialize();
    nd.N = 0;
    nd.temp_N = 0;
    nd.D = newD;
    nd.uniform = newUniform;
    nd.theta = newTheta;
    nd.Delta = newDelta;
    nd.offset = zeros(newD);
    nd.total = zeros(size(newTheta));
    nd.total2 = zeros(size(newTheta));
    nd.temp_total = zeros(size(newTheta));
    nd.temp_total2 = zeros(size(newTheta));
    nd.mean = zeros(size(newTheta));
    nd.sumlogDelta_div2 = sum(log(nd.Delta))/2;
end