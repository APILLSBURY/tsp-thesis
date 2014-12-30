function mean = NormalD_calc_mean(nd)
   if nd.uniform
       mean = nd.total ./ nd.N;
   else
       mean = (nd.theta + nd.offset + nd.Delta .* nd.total) ./ (1 + nd.N.*nd.Delta);
   end
end