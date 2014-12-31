function [] = compute_of(oim0, oim1, outname)

% add the optical flow path
addpath('flow_code_v2/');

if ~exist(outname, 'file')
    
    bv = estimate_flow_interface(oim1, oim0, 'classic+nl-fast');
    fv = estimate_flow_interface(oim0, oim1, 'classic+nl-fast');
    flow.bvx = bv(:, :, 1);
    flow.bvy = bv(:, :, 2);
    flow.fvx = fv(:, :, 1);
    flow.fvy = fv(:, :, 2);
    save(outname, 'flow');
end