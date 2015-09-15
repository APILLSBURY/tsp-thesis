% find the median theta in a mask
function [theta, theta_median, theta_offset] = get_theta_median(object_mask, vx, vy)
    if nargin==3
        theta = atan2(vy, vx);
        mask_theta = theta(object_mask);
    elseif nargin==1
        mask_theta = object_mask;
        theta = mask_theta;
    else
        error('Wrong number of arguments to get_theta_median');
    end
    curr_max = sum(mask_theta<0);
    quadrant = 1;
    if sum(mask_theta>-pi/2 & mask_theta<pi/2)>curr_max
        curr_max = sum(mask_theta>-pi/2 & mask_theta<pi/2);
        quadrant = 2;
    end
    if sum(mask_theta>0)>curr_max
        curr_max = sum(mask_theta>0);
        quadrant = 3;
    end
    if sum(mask_theta<-pi/2 | mask_theta>pi/2)>curr_max
        curr_max = sum(mask_theta<-pi/2 | mask_theta>pi/2);
        quadrant = 4;
    end
    
    if quadrant==1
        theta_offset = pi/2;
    elseif quadrant==2
        theta_offset = 0;
    elseif quadrant==3
        theta_offset = -pi/2;
    elseif quadrant==4
        theta_offset = pi;
    end
    theta = theta + theta_offset;
    theta(theta>pi) = theta(theta>pi)-2*pi;
    theta(theta<-pi) = theta(theta<-pi)+2*pi;
    if nargin==3
        theta_median = median(theta(object_mask));
    elseif nargin==1
        theta_median = median(theta);
    else
        error('Wrong number of arguments to get_theta_median');
    end
end