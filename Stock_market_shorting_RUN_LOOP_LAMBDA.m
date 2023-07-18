%Stock_market_shorting_RUN_LOOP_LAMBDA
%Stock market model with alternative uptick rule: policy analysis of kappa.
%This code uses loops over parameters to construct Figure 6 in the paper.
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk). This version: July 2023.

clc, clear, %close all

lambda_stack = [1 10000];   %First loss fig
kappa_star = NaN(length(lambda_stack),1);

for zz=length(lambda_stack):-1:1

    lambda = lambda_stack(zz);
        
    if zz == length(lambda_stack)
        Unconstrained = 1;
        run Stock_market_shorting_alt_uptick_IRF_LOOP
        Loss_unc = Loss; Loss_unc_1 = Loss_unc(1);
        x_abs_unc = x_abs(1); sum_Gini_unc = sum_Gini(1);
    end
        
        Unconstrained = 0;
        run Stock_market_shorting_alt_uptick_IRF_LOOP

        %if zz==length(lambda_stack)
        maxi = max(Loss);
        %end
        

[minL(zz),min_loc(zz)] = min(Loss);
kappa_star(zz) = kappa_stack(min_loc(zz));

    if minL(zz) == Loss(end)
        kappa_star(zz) = kappa_stack(end);
    end

       for jj = 1:length(kappa_stack)

       if jj==1.5 % jj===1|| jj==3 || jj==length(kappa_stack) || jj == round(length(kappa_stack)/2+1) %|| jj==min_loc(zz) 
       labels{jj} = round(kappa_stack(jj),4);
       else
          labels{jj} = [];
      end

   end

figure(1)
subplot(1,3,3), plot(kappa_stack,Loss/maxi, 'Marker', 'x', 'MarkerSize', 3, 'color', '[0.2 0.2 0.2]','LineWidth', 1.2),
title('Loss: L_\kappa'), axis([-inf,inf,-inf,inf]), xlabel('\kappa'), hold on,

    if zz == 1
        subplot(1,3,1), plot(kappa_stack,x_abs/x_abs_unc,'Marker', 'x', 'MarkerSize', 3, 'color', '[0.2 0.2 0.2]','LineWidth', 1.2), title('Mispricing'), 
        axis([-inf,inf,-inf,inf]), xlabel('\kappa'), hold on 
        subplot(1,3,2), plot(kappa_stack,sum_Gini/sum_Gini_unc, 'Marker', 'x', 'MarkerSize', 3, 'color', '[0.2 0.2 0.2]','LineWidth', 1.2), 
        axis([-inf,inf,-inf,inf]), xlabel('\kappa'), hold on, title('Inequality')
        
        figure(2)
        subplot(1,3,1), plot(sum_Gini/sum_Gini_unc,x_abs/x_abs_unc,'-o', 'color','[0.1 0.1 0.1]'), ylabel('Mispricing'), xlabel('Inequality') 
        hold on, text(sum_Gini/sum_Gini_unc,x_abs/x_abs_unc,labels,'VerticalAlignment','top','HorizontalAlignment','left'),
        title('\beta = 1.55'), %axis([-inf,inf,-inf,inf])
    end


end
