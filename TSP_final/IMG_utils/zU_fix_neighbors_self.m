function IMG_SP = U_fix_neighbors_self(IMG_SP, IMG_label, k)
    IMG_SP(k).neighbors = zeros(size(IMG_SP(k).neighbors));
    for i=find(IMG_SP(k).borders)'
        IMG_SP(k) = SP_update_neighbors_add_self(IMG_SP(k), IMG_label, i);
    end
end