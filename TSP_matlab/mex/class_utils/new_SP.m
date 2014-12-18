% --------------------------------------------------------------------------
% -- SP
% --   constructor; initializes to empty new super pixel
% --------------------------------------------------------------------------
function SP = new_SP(pos, app, UID, prev_v, max_N)
    SP.N = 0;
    SP.UID = UID;
    SP.pos = pos;
    SP.app = app;
    SP.prev_v = prev_v;
    SP.borders = false(max_N, 1);
    SP.pixels = false(max_N, 1);
    SP.neighbors = zeros(max_N, 1);
    SP.log_likelihood = NormalD_calc_logposterior(pos) + NormalD_calc_logposterior(app);
    SP.log_likelihood_empty = SP.log_likelihood;
end