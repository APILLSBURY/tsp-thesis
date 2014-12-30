function nd = NormalD_update_stats_rem(nd, data)
   nd.temp_N = nd.N - 1;
   nd.temp_total = nd.total - data;
   nd.temp_total2 = nd.total2 - data.^2;
end