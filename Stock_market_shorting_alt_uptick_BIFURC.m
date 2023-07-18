%Stock market model with short-selling constraint and endogenous shares: bifurcation diagrams 
%Last updated: July 26 2022. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk).

clear; clc; %close all; 

%------------------
%Parameter values
%------------------
H = 1000;  r = 0.1; a = 1; 
betta = 5; %betta = 0 gives case of fixed shares
dbar = 0.6; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r; %Fundamental price
kappa = 0.1;  %Alt. uptick rule

%----------------
%Coding choices
%----------------
Iter = 1;  %Iter = 1 turns on iterative algorithm (advisable for large H).
n_iter = 6; % no. of iterations 
Back = 0;  %Forward search is the default. Else solve backward from most optimistic type. 
Unconstrained = 1; %Set Unconstrained = 1 to simulate without short-selling constraints. 

%----------------------
%Specify belief types
%----------------------

run Stock_market_shorting_insert_benchmark
%run Stock_market_shorting_insert_disperse   %Old file to plot additional results
%run Stock_market_shorting_insert_disperse2   %Old file to plot additional results

num_betta = length(betta_stack);  dev = NaN(num_betta,1); dev1 = dev; dum = dev;

%-----------------
%Dividend shocks
%-----------------
shock = zeros(T,1);  %Deterministic skeleton

%------------------------
%Specify initial values
%------------------------
n_init = 1/H*ones(1,H);

%Baseline case
rng(3); init_stack = pf - 4*rand(M,1);

%Disperse beliefs
%rng(3);  init_stack = pf - 4*rand(M,1); 

%Disperse beliefs 2
%rng(3); init_stack = pf - 4 + 8*rand(M,1);

%----------------------
%Preallocate matrices
%----------------------
brk = zeros(M,1); percent = zeros(length(betta_stack),1); C1 = brk; C11 = brk; C12 = brk; sd_x = brk; 
x_stack = NaN(window,M); x_plot=NaN(window*M,1); Dem = NaN(H,T); sd_plot1=zeros(M,num_betta);
U = NaN(H,1); Bind = zeros(T,1);  D = NaN(H,1); Check1 = D; Check11 = D;

for v=1:num_betta 
    
    betta = betta_stack(v);
    x = NaN(T,1); 
    
for m = 1:M

%Initial price    
p0 = init_stack(m); x0 = p0 - pf; xlag = p0 - pf;
Bind = zeros(T,1);

for t=1:T 
    
    Beliefs = NaN(H,1);
    
    if t==1
        Beliefs = b + g*x0;
        n = n_init;
        cond = x0 - xlag + kappa*abs(xlag+pf); 
    elseif t==2
        Beliefs = b + g*x(t-1);
        n = n_init;
        cond = x(1) - x0 + kappa*abs(x0+pf); 
    elseif t>=3
        Beliefs = b + g*x(t-1);
        if t==3
            Dlag2 = (b + g*x0 + a*sigma^2*Zbar - (1+r)*x(t-2))/(a*sigma^2);
        else
            Dlag2 = (b + g*x(t-3) + a*sigma^2*Zbar - (1+r)*x(t-2))/(a*sigma^2);
        end
        if Bind(t-2) == 1
            Dlag2(Dlag2<0) = 0;
        end
        U = exp(betta*( (x(t-1) + a*sigma^2*Zbar + shock(t-1) - (1+r)*x(t-2))*Dlag2 - C) );
        n = transpose(U)/sum(U);
        cond = x(t-1) - x(t-2) + kappa*abs(x(t-2)+pf);
    end
    
        [Beliefs_sort,I] = sort(Beliefs);  
        n_adj = n(I);  %Sort pop. shares with beliefs
       
%Trial unconstrained solution
xstar = n*Beliefs/(1+r);   

if n*Beliefs - min(Beliefs) > a*sigma^2*Zbar && Unconstrained == 0 && cond <=0 
        
        Bind(t) = 1;
        
%Sort beliefs when there are ties (uncomment to use, not essential)
        %if length(unique(Beliefs)) ~= H
        %    Stock_market_shorting_sort_insert
        %end
            
%Obtain initial guess for no. short-sellers
        Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar)/(a*sigma^2);   
        k_init0 = sum(Demand_star<0);
   
%Decide whether to run iterations
   if Iter == 0
       k_init = k_init0;
   else 
       Stock_market_shorting_iterations_insert
       %num_iter(t) = n_iter;
   end
   
%------------------------------------------
%Find the equilibrium no.of short-sellers
%------------------------------------------
    %if Back == 0
        for k = k_init:length(Beliefs_sort)-1
            if n_adj(k:end)*Beliefs_sort(k:end) - sum(n_adj(k:end))*Beliefs_sort(k) > a*sigma^2*Zbar && n_adj(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj(k+1:end))*Beliefs_sort(k+1) <= a*sigma^2*Zbar 
            break
            end
        end
            
    %else    %uncomment for backward search
    %    for k = length(Beliefs_sort)-1:-1:k_init
     %       if n_adj(k:end)*Beliefs_sort(k:end) - sum(n_adj(k:end))*Beliefs_sort(k) > a*sigma^2*Zbar && n_adj(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj(k+1:end))*Beliefs_sort(k+1) <= a*sigma^2*Zbar 
     %       break
     %       end
     %  end
    %end 

kstar = k;  %Bind_no(t) = k;
          
       x(t) = ( n_adj(kstar+1:end)*Beliefs_sort(kstar+1:end) - sum(n_adj(1:kstar))*a*sigma^2*Zbar  ) / ( (1+r)*sum(n_adj(kstar+1:end)) );   
       
else 
        x(t) = xstar;    %Solution when SS constraints are slack or ignored       
end

%------------------------
%Check market clearing
%------------------------
    D = (Beliefs + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    D_adj = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    if Bind(t) == 1
        D(D<0) = 0;  
        D_adj(D_adj<0) = 0; 
    end
    Check1(t) = abs(n*D - Zbar); 
    Check11(t) = abs(n_adj*D_adj - Zbar);

end

%Store value for bifurc diagram
x_stack(1:window,m) = x(end+1-window:end); 
%Check for no attractor
r1 = 1-isreal(x(end)); r2 = isnan(x(end)); r3 = isinf(x(end));

%-------------------------------
%Record sims with no attractor
%-------------------------------
if (r1+r2+r3)>0
    brk(m) = 1;  
end

C1(m) = max(Check1);
C11(m) = max(Check11);
C12(m) = max(Bind);
sd_x(m) = std(x(end+1-50:end));

end

x_plot(:,v) = reshape(x_stack,1,[]);
sd_plot1(:,v) = sd_x;
%sd_plot1(:,v) = sd_plot1(:,v)./sd_plot(:,v);
%Uncomment after simulating Unconstrained = 0 to get relative SD for plot

%Percentage of sims with no attractor
percent(v) = 100*sum(brk)/M;

dev(v) = max(C1);
dev1(v) = max(C11);
dum(v) = max(C12);

end

%Accuracy checks
max(dev)
max(dev1)
%Check whether SS constraint binds in one or more sims
max(dum)

%-----------------
% Plot results
%-----------------
Bifurcation_plotter_base

%Bifurcation_plotter_comp  %Old file to plot additional results
%Bifurcation_plotter_comp2  %Old file to plot additional results
        
