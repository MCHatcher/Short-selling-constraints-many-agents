%Bifurcation_plotter_comp

%Parameter values
%H = 1000;  r = 0.1; a = 1; betta = 5; dbar = 0.6; sigma = 1; Zbar = 0.1;  
%pf = (dbar - a*sigma^2*Zbar)/r; kappa = 0.1; 

%Disperse beliefs
%M = 30;  %No. of initial values 
%window = 10;  %Sample to plot  
%T = 3000+window; rng(3); 
%num_betta = 82; betta_min = 0.9999; betta_max = 6.01; 

%b = zeros(H,1); C = b; g = 1.2*ones(H,1); 
%g(1:H/2) = 0; b(1:H/2) =  linspace(-0.2,0.2,H/2);  C(1:H/2) = 1-abs(b(1:H/2));

%betta_stack1 = linspace(betta_min,2.6,8); betta_stack2 = linspace(2.601,betta_max,num_betta);
%betta_stack = [betta_stack1 betta_stack2];
%num_betta = length(betta_stack);
%init_stack = pf - 4*rand(M,1); 

subplot(2,2,2)
hold on,
xlabel('Intensity of choice \beta'), ylabel('Price deviation \it{x}_t'), title('Alternative uptick rule'), %title('Alternative uptick rule')
axis([betta_min,betta_max,-2.3,10.3]), set(gca, 'box','on')

for v=1:num_betta
   
%plot(betta_stack(v), x_plot(:,v),'o', 'MarkerSize', 2.2, 'Color','k')  %[0.5,0.5,0.5]
plot(betta_stack(v), x_plot(:,v),'.', 'MarkerSize', 6.5, 'Color',[0.5,0.5,0.5])
%plot(betta_stack(1:20), max(x_plot(:,1:20)),'--','Linewidth', 4, 'Color',[1,1,1])

end

subplot(2,2,4), hold on

for v=1:num_betta

plot(betta_stack(v), sd_plot1(:,v),'.', 'MarkerSize', 6.5, 'Color','k')
%plot(betta_stack,sd_x,'Color',[0.5,0.5,0.5]), , 

end

title('Price volatility: zoomed in')
xlabel('Intensity of choice \beta'), ylabel('Volatility ratio'), axis([betta_min,betta_max,-inf,inf])

%Line(1:length(betta_stack)) = 0;
%plot(betta_stack,Line,'--k')

%subplot(2,2,3), plot(percent,'k')

