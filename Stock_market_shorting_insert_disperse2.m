%Stock_market_shorting_insert_disperse2

%Parameter values
M = 30;  %No. of initial values 
window = 15;  %Sample to plot  
T = 3000+window; rng(3);
num_betta = 82; betta_min = 0.999; betta_max = 6.01;

%Belief types
rng(10)
b = zeros(H,1); C = b; g = b; g(H/2+1:H) = 1 + 0.4*rand(H/2,1);
b(1:H/2) =  linspace(-0.2,0.2,H/2);  C(1:H/2) = 1-abs(b(1:H/2));

%Bifurcation paramater
betta_stack1 = linspace(betta_min,2.6,8); betta_stack2 = linspace(2.601,betta_max,num_betta);
betta_stack = [betta_stack1 betta_stack2];


