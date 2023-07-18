%Stock_market_shorting_RUN_LOOP_BETTA
%Stock market model with alternative uptick rule: policy analysis of kappa.
%This code uses loops over parameters to construct Figure 7 in the paper.
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). This version: July 2023.

clc, clear, %close all

betta_stack = linspace(3,4,13);  
lambda_stack = linspace(1,10000,2);  

kappa_star = NaN(length(lambda_stack),length(betta_stack));

for vv=1:length(lambda_stack)

for zz=1:length(betta_stack)

    betta = betta_stack(zz);
    lambda = lambda_stack(vv);
       
    Unconstrained = 0;
    run Stock_market_shorting_alt_uptick_IRF_LOOP

    [minL,min_loc] = min(Loss);
    kappa_star(vv,zz) = kappa_stack(min_loc);

    if minL == Loss(end)
        kappa_star(vv,zz) = kappa_stack(end);
    end

end

end

%---------------
%Plot results
%--------------

figure(1)
subplot(1,2,1), plot(betta_stack,kappa_star(1,:), 'Marker', 'x', 'MarkerSize', 4, 'color', '[0.2 0.2 0.2]','LineWidth', 1), axis([-inf,inf,0,0.1]),
title('\lambda = 1'), xlabel('Intensity of choice \beta'), ylabel('Optimal \kappa')
subplot(1,2,2), plot(betta_stack,kappa_star(2,:), 'Marker', 'x', 'MarkerSize', 4, 'color', '[0.2 0.2 0.2]','LineWidth', 1), axis([-inf,inf,0,0.1]),
title('\lambda = 10000'), xlabel('Intensity of choice \beta'), ylabel('Optimal \kappa')

