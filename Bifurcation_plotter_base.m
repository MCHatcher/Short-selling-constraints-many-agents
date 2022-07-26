%Bifurcation_plotter_base

%Disperse beliefs
%M = 20;  %No. of initial values 
%window = 15;  %Sample to plot  
%T = 3000+window; rng(3); 
%num_betta = 150; betta_min = 0.999; betta_max = 5.01; 

%Parameter values
%Two beliefs
%b = zeros(H,1); C = b; g = linspace(1.2,1.2,H); g = g'; g(1:H/2) = 0; C(g==0) = 1; 
%betta_stack1 = linspace(betta_min,3.85,20); betta_stack2 = linspace(3.851,betta_max,num_betta);
%betta_stack = [betta_stack1 betta_stack2];
%num_betta = length(betta_stack);

%Initial values
%init_stack = pf - 4*rand(M,1); 

figure(1)
%subplot(1,2,1)
hold on,
xlabel('Intensity of choice \beta'), ylabel('Price deviation \it{x}_t'), %title('Absence of short-selling constraints'), %title('Alternative uptick rule')
axis([betta_min,betta_max,-2.2,1]), set(gca, 'box','on')

for v=1:num_betta
   
plot(betta_stack(v), x_plot(:,v),'o', 'MarkerSize', 2.2, 'Color','k') %[0.5,0.5,0.5]
%plot(betta_stack(1:20), max(x_plot(:,1:20)),'--','Linewidth', 4, 'Color',[1,1,1])

end

%Line = -2.2 + 3.2*rand(500,1);
%coef(1:500) = 3.86;
%plot(coef,Line,'--k')


