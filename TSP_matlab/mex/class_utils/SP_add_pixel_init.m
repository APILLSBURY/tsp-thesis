% --------------------------------------------------------------------------
% -- add_pixel_init
% --   Adds a pixel to the super pixel, updates the appearance and position
% -- parameters, and the linked lists. Does not update the likelihoods.
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
function SP  = SP_add_pixel_init(SP, data, index, is_border, doApp)
    SP.N = SP.N+1;
    SP.pos = NormalD_add_data(SP.pos, data(index, 1:2));
    
    if (doApp)
        SP.app = NormalD_add_data(SP.app, data(index, 3:5));
    end
    
    SP.pixels(index) = true;
    
    if is_border
        SP.borders(index) = true;
    end
end