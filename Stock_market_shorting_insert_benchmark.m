%Stock_market_shorting_insert_benchmark. Used to plot Fig. 2 in the paper.

%Parameter values
M = 20;  %No. of initial values 
window = 15;  %Sample to plot  
T = 3000+window; rng(3); 
num_betta = 150; betta_min = 0.999; betta_max = 5.01;

%Two beliefs
b = zeros(H,1); C = b; g = linspace(1.2,1.2,H); g = g'; g(1:H/2) = 0; C(g==0) = 1; 

%Bifurcation parameter
betta_stack1 = linspace(betta_min,3.85,20); betta_stack2 = linspace(3.851,betta_max,num_betta);
betta_stack = [betta_stack1 betta_stack2];
num_betta = length(betta_stack);