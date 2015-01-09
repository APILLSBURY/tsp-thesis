function empty = SP_is_empty(IMG, index)
    empty = (index > numel(IMG_SP) || isempty(IMG_SP(index).N) || IMG_SP(index).N == 0);
end