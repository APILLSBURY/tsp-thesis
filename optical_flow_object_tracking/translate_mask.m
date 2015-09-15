function mask2 = translate_mask(mask, vx, vy)
    mask2 = false(size(mask));
    for r = 1:size(mask, 1)
        for c = 1:size(mask, 2)
            if mask(r, c)
                mask2(r + round(vy(r, c)), c + round(vx(r, c))) = true;
            end
        end
    end
end
            