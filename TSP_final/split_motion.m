vx = abs(flow.bvx);
vy = abs(flow.bvy);
stdx = mean(std(vx));
stdy = mean(std(vy));
meanx = mean(mean(vx));
meany = mean(mean(vy));
vx(vx<meanx+2*stdx) = 0;
vy(vy<meany+2*stdy) = 0;