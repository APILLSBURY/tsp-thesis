% --------------------------------------------------------------------------
% -- move_merge_calc_delta
% --   calculates the probability of merging the superpixel at k with the
% --   superpixel at merge_k.
% --
% --   parameters:
% --     - k : the original superpixel
% --     - merge_k : the superpixel to merge with
% --------------------------------------------------------------------------

function prob = U_move_merge_calc_delta(IMG, k, merge_k)
    prob = SP_log_likelihood_test_merge1(IMG.SP(k), IMG.SP(merge_k)) + IMG.SP(merge_k).log_likelihood_empty;
    prob = prob - IMG.SP(k).log_likelihood - IMG.SP(merge_k).log_likelihood;

    prob = prob + U_calc_model_order(IMG, IMG.SP(k).N + IMG.SP(merge_k).N, IMG.SP_old(k));
    prob = prob + U_calc_model_order(IMG, 0, IMG.SP_old(merge_k));
    prob = prob - U_calc_model_order(IMG, IMG.SP(k).N, IMG.SP_old(k));
    prob = prob - U_calc_model_order(IMG, IMG.SP(merge_k).N, IMG.SP_old(merge_k));
end