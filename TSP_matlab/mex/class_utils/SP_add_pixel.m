% --------------------------------------------------------------------------
% -- add_pixel
% --   Adds a pixel to the super pixel, updates the appearance and position
% -- parameters, the linked lists, and the likelihoods. Same as
% -- add_pixel_init except it also updates the likelihood.
% --
% --   parameters:
% --     - data : [5,N] matrix containing all the data in an image
% --     - index : in [0,N-1], index to a [5,1] data vector
% --     - is_border : indicator as to whether or not to add to border LL
% --     - doApp (true) : indicates if the appearance should be added also
% --   return parameters:
% --     - pixel_ptr : a pointer to the added linked list node pixel
% --     - border_ptr : a pointer to the added linked list node border
% --------------------------------------------------------------------------
function SP = SP_add_pixel(SP, data, index, is_border, doApp)
   SP = SP_add_pixel_init(SP, data, index, is_border, doApp);
   SP = SP_calculate_log_probs(SP);
end