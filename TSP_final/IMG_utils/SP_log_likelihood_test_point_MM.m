function prob = SP_log_likelihood_test_point_MM(SP, data)
    prob = NormalD_calc_logposterior_MM(SP.pos, data(:, 1:2)) + NormalD_calc_logposterior_MM(SP.app, data(:, 3:5));
end