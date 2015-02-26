% =============================================================================
% == merge_move.m
% == --------------------------------------------------------------------------
% == Finds the optimal merge between all pairs of super pixels.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-4-2014
% =============================================================================

function [IMG_K, IMG_label, IMG_SP, IMG_SP_changed, IMG_alive_dead_changed, IMG_SP_old]  = merge_move(IMG_label, IMG_SP, IMG_SP_old, IMG_alive_dead_changed, IMG_SP_changed, model_order_params, IMG_K, its)
    disp('merge_move');
    for i=1:its
        % choose a random order of super pixels
        perm = randperm(IMG_K);


        for kindex=1:length(perm)
            k=perm(kindex);

            % find a nonempty super pixel
            if ~(k > numel(IMG_SP) || isempty(IMG_SP(k).N) || IMG_SP(k).N == 0)
                [IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed] = merge_move_internal(IMG_label, IMG_SP, IMG_K, IMG_SP_changed, IMG_alive_dead_changed, IMG_SP_old, model_order_params, k, false);
            end
        end
    end
end



