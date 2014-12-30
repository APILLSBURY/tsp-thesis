function nd = NormalD_merge(nd, nd2)
    nd.N = nd.N + nd2.N;
    nd.total = nd.total + nd2.total;
    nd.total2 = nd.total2 + nd2.total2;
end