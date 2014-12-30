function nd = NormalD_update_stats_new2(nd, nd_other1, nd_other2)
   nd.temp_N = nd.N + nd_other1.N + nd_other2.N;
   nd.temp_total = nd.total + nd_other1.total + nd_other2.total;
   nd.temp_total2 = nd.total2 + nd_other1.total2 + nd_other2.total2;
end