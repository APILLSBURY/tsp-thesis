function empty = SP_is_empty(IMG, index)
    empty = index > numel(IMG.SP) || isempty(IMG.SP(index).N) || IMG.SP(index).N == 0;
end