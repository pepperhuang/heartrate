function bayes_prc_fit(my_id,~,to_dir)

% Load PRC data
load([to_dir 'phase_data/' my_id '_prc.mat'],'prc_freq','prc_scramble');

% Function to calculate the sum of squared errors between the current
% parametrized fit and the observed PRC
function res = fun(x,y)
    res = ones(1,24)*x(4);
    xvals = mod((1:24)-(x(1)-0.5*x(2)),24);
    xvals_L = (xvals<x(2));
    xvals_R = ~xvals_L;
    res(xvals_L) = res(xvals_L) + x(3)*sin((pi/x(2))*(xvals(xvals_L)-0.5*x(2)));
    res(xvals_R) = res(xvals_R) - x(3)*sin((pi/(24-x(2)))*(xvals(xvals_R)-12-0.5*x(2)));
    res = sum((res-y').^2);
end
fun_handle = @(x) fun(x,-prc_freq);

% Fit the parametrized curve by minimizing the sum of squared errors
myfit = fminsearchbnd(fun_handle,[12 6 2 0],[-inf 0 0 -inf],[inf 24 inf inf]);

% The first parameter (phase shift) is periodic
myfit(1) = mod(myfit(1),24);

% Plot the actual PRC data
scatter(-11:12,-prc_freq)

% Compare with the parametrized fit
hold on;
res = ones(1,24)*myfit(4);
nxvals = mod((1:24)-(myfit(1)-0.5*myfit(2)),24);
nxvals_L = (nxvals<myfit(2));
nxvals_R = ~nxvals_L;
res(nxvals_L) = res(nxvals_L) + myfit(3)*sin((pi/myfit(2))*(nxvals(nxvals_L)-0.5*myfit(2)));
res(nxvals_R) = res(nxvals_R) - myfit(3)*sin((pi/(24-myfit(2)))*(nxvals(nxvals_R)-12-0.5*myfit(2)));
plot(-11:12,res);
hold off;

% Title and save plot
title([num2str(myfit(1)) ', ' num2str(myfit(2)) ', ' num2str(myfit(3)) ', ' num2str(myfit(4))])
save([to_dir 'phase_data/' my_id '_prc.mat'],'prc_freq','prc_scramble','myfit');
end