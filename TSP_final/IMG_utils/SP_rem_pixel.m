% --------------------------------------------------------------------------
% -- rem_pixel
% --   Removes a pixel from the super pixel, updates the appearance and
% -- position parameters, linked lists, and likelihoods.
% --
% --   parameters:
% --     - data : [5,N] matrix containing all the data in an image
% --     - index : in [0,N-1], index to a [5,1] data vector
% --     - pixel_ptr : a poitner into the pixels linked list to be removed
% --     - border_ptr : a pointer into the borders linked list to be removed
% --     - doApp (true) : indicates if the appearance should be done also
% --------------------------------------------------------------------------
function SP = SP_rem_pixel(SP, data, index, doApp)
   if (SP.N<=0)
      disp('Trying to remove a pixel from an empty super pixel!');
   end

   SP.N = SP.N-1;
   SP.pos = NormalD_rem_data(SP.pos, data(index, 1:2));
   
   if (doApp)
      SP.app = NormalD_rem_data(SP.app, data(index, 3:5));
   end

   SP.borders(index) = false;
   SP.pixels(index) = false;

   SP = SP_calculate_log_probs(SP);
end
