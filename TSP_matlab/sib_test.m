num_2_circles = 0;
num_3_circles = 0;
max_iters = 10000;
for iter=1:max_iters
    %find an order where nobody gives to themself
    has_not_given = true(1, 7);
    giving_to = randperm(7);
    while any(giving_to==1:7)
        giving_to = randperm(7);
    end
    
    %check to see if there are two circles
    sib = 1;
    while has_not_given(sib)
        has_not_given(sib) = false;
        sib = giving_to(sib);
    end
    if any(has_not_given)
        num_2_circles = num_2_circles + 1;
    end
    
    %check for a third circle
    sib = find(has_not_given, 1);
    while has_not_given(sib)
        has_not_given(sib) = false;
        sib = giving_to(sib);
    end
    if any(has_not_given)
        num_3_circles = num_3_circles + 1;
    end
end

disp(num_2_circles/max_iters);
disp(num_3_circles/max_iters);
        
