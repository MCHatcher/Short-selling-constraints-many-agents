%Stock market model with short-selling constraint and endogenous shares: simulations 
%Simulations of price scenarios (Fig. 3), omitting simulation of wealth
%This code is faster than the ...SIMS_fast.m version when the number of
%types is very large because belief dispersion is computed recursively
%Feb 6 2024. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clear, clc, %close all; 

%------------------
%Parameter values
%------------------
H = 1000000;  %No. of types
r = 0.1; a = 1; 
betta = 4.5; %3, 4.5
dbar = 0.6; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r; %Fundamental price
kappa = 0;  %Alternative uptick rule

%----------------
%Coding choices
%----------------
Iter = 1;  %Iter = 1 turns on iterative algorithm (advisable for large H).
Fixed = 0; %Fixed  = 1: Pick fixed rather than time-varying (fitness-based) population shares. 
n_iter = 6; % no. of iterations (increase for large H)
Unconstrained = 0; %Set Unconstrained = 1 to simulate without short-selling constraints. 
T = 500;  %no. of periods

%----------------------
%Preallocate matrices
%----------------------
U = NaN(H,1); Bind = zeros(T,1); Bind_no = Bind;  D = NaN(H,1); D_lag2 = D; x = NaN(T,1); Time = x;
Check1 = NaN(T,1); Check11 = Check1; Beliefs = NaN(H,1);

%--------------------------
%Generate dividend shocks 
%--------------------------
%Uncomment initially to store shocks in memory
rng(1), sigma_d  = 0.0099;   
pd = makedist('Normal','mu',0,'sigma',sigma_d);  %Truncated normal distribution
pd_t = truncate(pd,-dbar,dbar);
rng(1), shock = random(pd_t,T,1);   
%NB. For timed sims, store shock in workspace and comment out other lines

%-------------------------------
%Initial values and predictors 
%-------------------------------
p0 = 8; x0 = p0 - pf; xlag = p0 - pf; 
n_init = 1/H*ones(1,H); 

%Disperse beliefs (Scenario 3)
b = zeros(H,1); C = b; g = 1.2*ones(H,1); 
g(1:H/2) = 0; b(1:H/2) =  linspace(-0.2,0.2,H/2);  C(1:H/2) = 1-abs(b(1:H/2));

for t=1:T 
    
    %Beliefs = NaN(H,1);
    
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
       
%------------------------------    
%Trial unconstrained solution
%------------------------------
xstar = n*Beliefs/(1+r);   

        [Beliefs_sort,I] = sort(Beliefs);  
        n_adj = n(I);

if n*Beliefs - min(Beliefs) > a*sigma^2*Zbar && Unconstrained == 0 && cond <=0 
        
        Bind(t) = 1;
        
%Sort beliefs when there are ties (uncomment to use, not essential)
        %if length(unique(Beliefs)) ~= H
        %    run Stock_market_shorting_sort_insert
        %end
    
%--------------------------------------------       
%Obtain initial guess for no. short-sellers
%--------------------------------------------
   Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar)/(a*sigma^2);    
   k_init0 = sum(Demand_star<0);
   
   Stock_market_shorting_iterations_insert 
   
%-----------------------------------------   
%Find the equilibrium no.of short-sellers
%-----------------------------------------
   
   disp = n_adj(k_init:end)*Beliefs_sort(k_init:end) - sum(n_adj(k_init:end))*Beliefs_sort(k_init);
   sum_n = sum(n_adj(k_init:end));

       for k = k_init:length(Beliefs_sort)-1
           sum_n = sum_n - n_adj(k);
           disp_init = disp;
           disp = disp - sum_n*(Beliefs_sort(k+1)- Beliefs_sort(k));

            if disp <= a*sigma^2*Zbar && disp_init > a*sigma^2*Zbar
                break
            end
       end

        %Original approach
        %for k = k_init:length(Beliefs_sort)-1
        %    if n_adj(k:end)*Beliefs_sort(k:end) - sum(n_adj(k:end))*Beliefs_sort(k) > a*sigma^2*Zbar && n_adj(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj(k+1:end))*Beliefs_sort(k+1) <= a*sigma^2*Zbar 
        %    break
        %    end
        %end

kstar = k;  Bind_no(t) = k;   %No. of constrained types
          
       x(t) = ( n_adj(kstar+1:end)*Beliefs_sort(kstar+1:end) - sum(n_adj(1:kstar))*a*sigma^2*Zbar  ) / ( (1+r)*sum(n_adj(kstar+1:end)) );   
       
else 
        x(t) = xstar;   %Solution when SS constraints are slack or ignored           
end

%-------------------------------------------------------
%Check market clearing (uncomment for acccuracy checks)
%-------------------------------------------------------
    D = (Beliefs + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    D_adj = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    if Bind(t) == 1
        D(D<0) = 0;  
        D_adj(D_adj<0) = 0; 
    end
    Check1(t) = abs(n*D - Zbar); 
    Check11(t) = abs(n_adj*D_adj - Zbar);

end

%Accuracy checks
max(Check1)
max(Check11)
     


