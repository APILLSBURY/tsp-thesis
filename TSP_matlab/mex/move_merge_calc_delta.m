% --------------------------------------------------------------------------
% -- move_merge_calc_delta
% --   calculates the probability of assigning the pixel at index to the
% -- cluster of nindex (neighbor index).
% --
% --   parameters:
% --     - index : the new point to add
% --     - nindex : the neighbor to add this point to
% --------------------------------------------------------------------------

function prob = move_merge_calc_delta(IMG, k, merge_k)
    addpath('mex/class_utils/');
    prob = SP_log_likelihood_test_merge1(IMG.SP(k), IMG.SP(merge_k)) + IMG.SP(merge_k).log_likelihood_empty;
    prob = prob - IMG.SP(k).log_likelihood + IMG.SP(merge_k).log_likelihood;

    prob = prob + U_calc_model_order(IMG, IMG.SP(k).N + IMG.SP(merge_k).N, IMG.SP_old(k));
    prob = prob + U_calc_model_order(IMG, 0, IMG.SP_old(merge_k));
    prob = prob - U_calc_model_order(IMG, IMG.SP(k).N, IMG.SP_old(k));
    prob = prob - U_calc_model_order(IMG, IMG.SP(merge_k).N, IMG.SP_old(merge_k));
end