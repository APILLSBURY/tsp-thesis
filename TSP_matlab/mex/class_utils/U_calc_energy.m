function E = U_calc_energy(IMG)
    E = 0;
    numEmpty = 0;
    for k=1:IMG.N
        if ~SP_is_empty(IMG, k)
            IMG.SP(k) = SP_calculate_log_probs(IMG.SP(k));
            E = E + IMG.SP(k).log_likelihood + U_calc_model_order(IMG, IMG.SP(k).N, IMG.SP_old(k));
        else
            numEmpty = numEmpty + 1;
        end
    end
    new_SP = create_SP_new(IMG);
    E = E + numEmpty*(new_SP.log_likelihood + U_calc_model_order(IMG, 0, false));
end