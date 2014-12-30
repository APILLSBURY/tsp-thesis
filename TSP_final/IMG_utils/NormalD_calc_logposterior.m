% --------------------------------------------------------------------------
% -- calc_logposterior_internal
% --   calculate log p(mu|theta)p(x|mu) with either the optimal iid mu, or
% -- the current value of mu depending on fixedMean
% --------------------------------------------------------------------------
function logprob = NormalD_calc_logposterior(nd, useTemp)
    if nargin==1
        useTemp = false;
    end
    if useTemp
        this_total = nd.temp_total;
        this_total2 = nd.temp_total2;
        this_N = nd.temp_N;
    else
        this_total = nd.total;
        this_total2 = nd.total2;
        this_N = nd.N;
    end
    
    logprob = -nd.D*(this_N+0.5)*0.5*1.837877066409 - nd.sumlogDelta_div2;
    if (this_N>0)
        if nd.uniform
            logprob = logprob + sum(this_total.^2 ./ this_N - this_total2)/2;
        else
            logprob = logprob + sum((nd.Delta.*(this_total.^2) + 2*this_total.*(nd.theta+nd.offset) ... 
                - this_N*((nd.theta+nd.offset).^2)) / (2*this_N.*nd.Delta+2) - this_total2./2);
        end
    end
end