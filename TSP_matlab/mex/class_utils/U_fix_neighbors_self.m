function IMG = U_fix_neighbors_self(IMG, k)
    IMG.SP(k).neighbors = zeros(size(IMG.SP(k).neighbors));
    for i=find(IMG.SP(k).borders)'
        IMG.SP(k) = SP_update_neighbors_add_self(IMG.SP(k), IMG.label, i);
    end
end