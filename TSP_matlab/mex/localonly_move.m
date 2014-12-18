% =============================================================================
% == localonly_move.m
% == --------------------------------------------------------------------------
% == An interface to perform local TSP moves without flow estimation.
% == See m files for calling convention.
% ==
% == All work using this code should cite:
% == J. Chang, D. Wei, and J. W. Fisher III. A Video Representation Using
% ==    Temporal Superpixels. CVPR 2013.
% == --------------------------------------------------------------------------
% == Written in C++ by Jason Chang and Donglai Wei 06-20-2013
% == Converted to MATLAB by Andrew Pillsbury 12-4-2014
% =============================================================================

function IMG = localonly_move(IMG, its)
    addpath('mex/class_utils/');
    for i=1:its
        fprintf('its=%d\n', i);
        fprintf('making moves mod\n');
        IMG = local_move_internal(IMG);
        if ~IMG.changed
            break;
        end
    end
end
