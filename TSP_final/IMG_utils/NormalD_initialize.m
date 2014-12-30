% --------------------------------------------------------------------------
% -- NormalD
% --   constructor; initializes to empty
% --------------------------------------------------------------------------
function nd = NormalD_initialize()
    nd.D = 0;
    nd.offset = [];
    nd.theta = [];
    nd.Delta = [];
    nd.total = [];
    nd.total2 = [];
    nd.N = 0;
    nd.temp_total = [];
    nd.temp_total2 = [];
    nd.temp_N = 0;
    nd.mean = [];
    nd.uniform = false;
    nd.sumlog_Delta_div2 = 0;
end