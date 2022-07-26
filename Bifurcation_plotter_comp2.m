%Bifurcation_plotter_comp2

%Disperse beliefs 2
%M = 30;  %No. of initial values 
%window = 15;  %Sample to plot  
%T = 3000+window; rng(3); 
%num_betta = 82; betta_min = 0.999; betta_max = 6.01; 

%rng(10)
%b = zeros(H,1); C = b; g = b; g(H/2+1:H) = 1 + 0.4*rand(H/2,1);
%b(1:H/2) =  linspace(-0.2,0.2,H/2);  C(1:H/2) = 1-abs(b(1:H/2));

%betta_stack1 = linspace(betta_min,2.6,8); betta_stack2 = linspace(2.601,betta_max,num_betta);
%betta_stack = [betta_stack1 betta_stack2];
%betta_stack = linspace(betta_min, betta_max, num_betta);
%num_betta = length(betta_stack);

%rng(3);
%init_stack = pf - 4 + 8*rand(M,1);

subplot(2,2,2)
hold on,
xlabel('Intensity of choice \beta'), ylabel('Price deviation \it{x}_t'), title('Alternative uptick rule'), %title('Alternative uptick rule')
axis([betta_min,betta_max,-2,2.2]), set(gca, 'box','on')

for v=1:num_betta
   
%plot(betta_stack(v), x_plot(:,v),'o', 'MarkerSize', 2.2, 'Color',[0.5,0.5,0.5])  %[0.5,0.5,0.5]
plot(betta_stack(v), x_plot(:,v),'.', 'MarkerSize', 6.5, 'Color',[0.5,0.5,0.5])
%plot(betta_stack(1:20), max(x_plot(:,1:20)),'--','Linewidth', 4, 'Color',[1,1,1])

end

subplot(2,2,3), hold on, plot(betta_stack,percent,'Color',[0.5,0.5,0.5]), title('Price instability: percent of simulations'), 
xlabel('Intensity of choice \beta'), ylabel('% explosive / no attractor'), axis([betta_min,betta_max,-inf,inf])

%Line(1:length(betta_stack)) = 0;
%plot(betta_stack,Line,'--k')

%subplot(2,2,3), plot(percent,'k')

