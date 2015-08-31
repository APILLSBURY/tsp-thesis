% find the absolute theta given a relative theta and an offset
function theta = get_absolute_theta(theta, theta_offset)
    theta = theta - theta_offset;
    if theta>pi
        theta = theta - 2*pi;
    end
    if theta<-pi
        theta = theta + 2*pi;
    end
end