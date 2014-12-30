% --------------------------------------------------------------------------
% -- empty
% --   Empties out the SP
% --------------------------------------------------------------------------
function SP = SP_empty(SP)
    NormalD_empty(SP.pos);
    NormalD_empty(SP.app);

    SP.N = 0;
    SP.borders = false(size(SP.borders));
    SP.pixels = false(size(SP.pixels));
    SP.neighbors = zeros(size(SP.neighbors));
    SP = SP_calculate_log_probs(SP);
end