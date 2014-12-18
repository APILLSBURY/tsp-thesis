function IMG = U_fix_neighbors_neighbors(IMG, k, kignore)
    if narg_in < 3
        kignore = 0;
    end
    for k2=find(IMG.SP(k).neighbors)'
        if k2 ~= kignore
            IMG = U_fix_neighbors_self(IMG, k2);
        end
    end
end