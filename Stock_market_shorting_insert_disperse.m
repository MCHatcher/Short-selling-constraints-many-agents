%Stock_market_shorting_insert_disperse

%Parameter values
M = 30;  %No. of initial values 
window = 10;  %Sample to plot  
T = 3000+window; rng(3); 
num_betta = 82; betta_min = 0.9999; betta_max = 6.01; 
%betta_stack = linspace(betta_min, betta_max, num_betta);

%Disperse beliefs
b = zeros(H,1); C = b; g = 1.2*ones(H,1); 
g(1:H/2) = 0; b(1:H/2) =  linspace(-0.2,0.2,H/2);  C(1:H/2) = 1-abs(b(1:H/2));

%Bifurcation parameter
betta_stack1 = linspace(betta_min,2.6,8); betta_stack2 = linspace(2.601,betta_max,num_betta);
betta_stack = [betta_stack1 betta_stack2];