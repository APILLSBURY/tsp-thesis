function prob = SP_log_likelihood_test_merge2(SP, SP_other1, SP_other2)
    prob = NormalD_calc_logposterior_new2(SP.pos, SP_other1.pos, SP_other2.pos) + NormalD_calc_logposterior_new2(SP.app, SP_other1.app, SP_other2.app);
end