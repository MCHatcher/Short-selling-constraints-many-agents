%Bifurcation_plotter_base
%Plots the bifurcation diagram in Figure 2

figure(1)
%subplot(1,2,1)
hold on,
xlabel('Intensity of choice \beta'), ylabel('Price deviation \it{x}_t'), %title('Absence of short-selling constraints'), %title('Alternative uptick rule')
axis([betta_min,betta_max,-2.2,1]), set(gca, 'box','on')

for v=1:num_betta
   
plot(betta_stack(v), x_plot(:,v),'o', 'MarkerSize', 2.2, 'Color','k') %[0.5,0.5,0.5]
%plot(betta_stack(1:20), max(x_plot(:,1:20)),'--','Linewidth', 4, 'Color',[1,1,1])

end



