x=(1:10000000)*2;
tic;
for index=x
    temp=index;
end
toc;
tic;
for index=1:length(x)
    temp=x(index);
end
toc;