% =============================================================================
% == calc_sdf.m
% == --------------------------------------------------------------------------
% == The MEX interface file to calculate a signed distance function
% == --------------------------------------------------------------------------
% == Copyright 2011. MIT. All Rights Reserved.
% == Written by Jason Chang 06-13-2011 in C++
% == Converted to MATLAB by Andrew Pillsbury 12-3-2014
% =============================================================================

% --------------------------------------------------------------------------
% -- calc_sdf
% --   calculates a signed distance function for distanceIn (which only must
% -- contain the correct sign). Pixel i is only updated if maskIn[i]=true,
% -- and in a narrow band region if size band_size.
% --
% --   parameters:
% --     - (xdimDouble, ydimDouble) : the dimensions of the image
% --     - distanceIn : the new level set function to calculate the SDF of.
% --         Only the signs of this matter.
% --     - maskIn : the boolean array indicating which pixels to update
% --     - oldDistance : the old values of the SDF
% --     - band_size : the size of the narrow band
% --
% --   return parameters:
% --     - distance : the calculated SDF
% --------------------------------------------------------------------------
function distance = calc_sdf(distanceIn, maskIn, oldDistance, band_size, distance)
   [xdim, ydim] = size(oldDistance);

   unmarked = NaN(xdim, ydim);

   accepted = false(xdim, ydim);

   for x=1:xdim
      for y=1:ydim
          if (~maskIn(x, y))
             accepted(x, y) = true;
             distance(x, y) = oldDistance(x, y);
          elseif (((distanceIn(x, y) >= 0) && ((x>1 && distanceIn(x-1, y)<0) || ...
                                               (x<xdim && distanceIn(x+1, y)<0) || ...
                                               (y>1 && distanceIn(x, y-1)<0) || ...
                                               (y<ydim && distanceIn(x, y+1)<0))) || ...
                   ((distanceIn(x, y) < 0 ) && ((x>1 && distanceIn(x-1, y)>=0) || ...
                                               (x<xdim && distanceIn(x+1, y)>=0) || ...
                                               (y>1 && distanceIn(x, y-1)>=0) || ...
                                               (y<ydim && distanceIn(x, y+1)>=0))))

             accepted(x, y) = true;
             distance(x, y) = 0;
             unmarked(x, y) = 0.5;
             if (distanceIn(x, y) < 0)
                unmarked(x, y) = -0.5;
             end
          else
             if (distanceIn(x, y)<0)
                distance(x, y) = -band_size;
             else
                distance(x, y) = band_size;
             end
             accepted(x, y) = false;
          end
      end
   end

   unmarked_sorted =zeros(xdim*ydim, 3);
   for x=1:xdim
       for y=1:ydim
           index = (x-1)*ydim + y;
           unmarked_sorted(index,:) = [unmarked(x, y) x y];
       end
   end
   unmarked_sorted = sortrows(unmarked_sorted);
   % recursively pick off the smallest distance until all pixels are done
   i = 1;
   while (~isnan(unmarked_sorted(i, 1)) && i <= size(unmarked_sorted, 1))
      x = unmarked_sorted(i, 2);
      y = unmarked_sorted(i, 3);
      distance(x, y) = unmarked_sorted(i, 1);
      accepted(x, y) = true;

      if (x<xdim && ~accepted(x+1, y) )
          [distance, unmarked] = calc_dist(x+1, y, unmarked, accepted, band_size, distance);
      end
      if (x-1>0 && ~accepted(x-1,y) )
          [distance, unmarked] = calc_dist(x-1, y, unmarked, accepted, band_size, distance);
      end
      if (y<ydim && ~accepted(x, y+1) ) 
          [distance, unmarked] = calc_dist(x, y+1, unmarked, accepted, band_size, distance);
      end
      if (y-1>=0 && ~accepted(x, y-1) )
          [distance, unmarked] = calc_dist(x, y-1, unmarked, accepted, band_size, distance);
      end
   end
end

% --------------------------------------------------------------------------
% -- calc_dist
% --   calculates distance to zero level set at far point (x,y)
% --
% --   parameters:
% --     - (x,y) : the coordinate where the distance is needed
% --     - (xNew, yNew) : the newly accepted distance that touches (x,y)
% --     - unmarked : the tree of points that have a finite distance but
% --         have not been accepted
% --     - accepted : the 2D boolean matrix saying if the pixel is accepted
% --     - tree_ptr : a 2D matrix of pointers pointing to the node in the
% --         unmarked tree (if exists) or NULL (if doesn't exist)
% --     - band_size : the narrow band size
% --
% --   return parameters:
% --     - distance : the SDF, updated (x,y)
% --------------------------------------------------------------------------
function [distance, unmarked] = calc_dist(x, y, unmarked, accepted, band_size, distance)

   [xdim, ydim] = size(accepted);
   % do x direction
   a = 0;
   b = 0;
   c = -1;
   
   checkbck = (x-1>=0) && accepted(x-1, y);
   checkfwd = (x+1<xdim) && accepted(x+1, y);
   dist = 10*band_size;
   if (checkbck && ~checkfwd)
      dist = distance(x-1, y);
   elseif (~checkbck && checkfwd)
      dist = distance(x+1, y);
   elseif (checkbck && checkfwd)
      if (distance(x+1, y) > distance(x-1, y))
          dist = distance(x-1, y);
      else
          dist = distance(x+1, y);
      end
   end
   if (dist ~= 10*band_size)
      a = a + 1;
      b = b - 2*dist;
      c = c + dist*dist;
   end

   % do y direction
   checkbck = (y-1>=0) && accepted(x, y-1);
   checkfwd = (y+1<ydim) && accepted(x, y+1);
   dist = 10*band_size;
   if (checkbck && ~checkfwd)
      dist = distance(x, y-1);
   elseif (~checkbck && checkfwd)
      dist = distance(x, y+1);
   elseif (checkbck && checkfwd)
      if (distance(x, y+1) > distance(x, y-1))
          dist = distance(x, y-1);
      else
          dist = distance(x, y+1);
      end
   end
   if (dist ~= 10*band_size)
      a = a + 1;
      b = b - 2*dist;
      c = c + dist*dist;
   end

   
   if distance(x, y) > 0
       sign = 1;
   else
       sign = 0;
   end
   dist = (-b + sign * sqrt(b*b - 4*a*c)) / (2*a);

   if (abs(dist) < abs(distance(x, y)) && abs(dist)<=band_size)
      % add it to the tree because it had never been added before
      unmarked(x, y) = dist;
      distance(x, y) = dist;
   end
end

