function SP = SP_calculate_log_probs(SP)
    if SP.N<=0
        SP.log_likelihood = SP.log_likelihood_empty;
    else
        SP.log_likelihood = NormalD_calc_logposterior(SP.pos) + NormalD_calc_logposterior(SP.app);
    end
end