%Stock market model with short-selling constraint and fixed pop. shares: simulations 
%Last updated: Dec 12, 2022. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clear, clc, %close all; 

%------------------
%Parameter values
%------------------
H = 3000;
r = 0.1; a = 1; 
dbar = 1.1; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r; %Fundamental price
cbar = 1;
r = r-cbar; %Adjustment for demand function

%----------------
%Coding choices
%----------------
Iter = 1;  %Iter = 1 turns on iterative algorithm (advisable for large H). 
n_iter = 6; % no. of iterations (increase for large H)
Unconstrained = 0; %Set Unconstrained = 1 to simulate without short-selling constraints. 
T = 500;  %no. of periods
rng(5)  %Set random seed

%----------------------
%Preallocate matrices
%----------------------
Bind = zeros(T,1); Bind_no = Bind;  D = NaN(H,1); D_lag2 = D; p = NaN(T,1); Time = p;
Check1 = NaN(T,1); Check11 = Check1; Beliefs = NaN(H,1);

%-------------------------------
%Initial values and predictors 
%-------------------------------
p0 = pf + 0.6; plag = p0; plag1 = p0; plag2 = p0;
n_init = 1/H*ones(1,H); 
n = n_init;

%Disperse beliefs
b = zeros(H,1);  g1 = 0.5*rand(H,1); g2 = 0.2*rand(H,1); 
g3 = -0.1*rand(H,1); g4 = -0.1*rand(H,1);
gf = -0.2 -0.6*rand(H,1); u_stack = 0.04*randn(H,T); u_stack(1:H,1:10) = 0;

%Trend followers
g3(1:H/3) = 0; g4(1:H/3) = 0; gf(1:H/3) = 0; %u_stack(1:H/3,1:T) = 0;

%Contrarians
g1(H/3+1:H*2/3) = 0; g2(H/3+1:H*2/3) = 0; gf(H/3+1:H*2/3) = 0; %u_stack(H/3+1:H*2/3,1:T) = 0;

%Fundamentalists (Arbitrageurs)
g1(H*2/3+1:end) = 0; g2(H*2/3+1:end) = 0; g3(H*2/3+1:end) = 0; g4(H*2/3+1:end) = 0;

for t=1:T 
    
    %Beliefs = NaN(H,1);
    
    u = u_stack(:,t);
    
    if t==1
        Beliefs = b + (g1+g3)*(p0 - plag) + (g2+g4)*(plag-plag2) + gf*(p0-pf) + u + dbar -a*sigma^2*Zbar; 
    elseif t==2
        Beliefs = b + (g1+g3)*(p(t-1) - p0) + (g2+g4)*(p0 -plag) + gf*(p(t-1)-pf) + u + dbar -a*sigma^2*Zbar;
    elseif t==3
        Beliefs = b + (g1+g3)*(p(t-1) - p(t-2)) + (g2+g4)*(p(t-2) -p0) + gf*(p(t-1)-pf) + u + dbar -a*sigma^2*Zbar;
    else
        Beliefs = b + (g1+g3)*(p(t-1) - p(t-2)) + (g2+g4)*(p(t-2) -p(t-3)) + gf*(p(t-1)-pf) + u + dbar -a*sigma^2*Zbar;
    end
 
       
%------------------------------    
%Trial unconstrained solution
%------------------------------
pstar = n*Beliefs/(1+r);   

        [Beliefs_sort,I] = sort(Beliefs);  
        n_adj = n(I);

if n*Beliefs - min(Beliefs) > a*sigma^2*Zbar && Unconstrained == 0 
        
        Bind(t) = 1;
        
%Sort beliefs when there are ties (uncomment to use, not essential)
        %if length(unique(Beliefs)) ~= H
        %    run Stock_market_shorting_sort_insert
        %end
    
%--------------------------------------------       
%Obtain initial guess for no. short-sellers
%--------------------------------------------
   Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*pstar)/(a*sigma^2);    
   k_init0 = sum(Demand_star<0);
   
   %run Stock_market_shorting_iterations_insert
   Stock_market_shorting_iterations_insert 
   
%-----------------------------------------   
%Find the equilibrium no.of short-sellers
%-----------------------------------------
        for k = k_init:length(Beliefs_sort)-1
            if n_adj(k:end)*Beliefs_sort(k:end) - sum(n_adj(k:end))*Beliefs_sort(k) > a*sigma^2*Zbar && n_adj(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj(k+1:end))*Beliefs_sort(k+1) <= a*sigma^2*Zbar 
            break
            end
        end

kstar = k;  Bind_no(t) = k;   %No. of constrained types
          
       p(t) = ( n_adj(kstar+1:end)*Beliefs_sort(kstar+1:end) - sum(n_adj(1:kstar))*a*sigma^2*Zbar  ) / ( (1+r)*sum(n_adj(kstar+1:end)) );   
       
else 
        p(t) = pstar;   %Solution when SS constraints are slack or ignored          
end

%-------------------------------------------------------
%Check market clearing (uncomment for acccuracy checks)
%-------------------------------------------------------
    D = (Beliefs + a*sigma^2*Zbar - (1+r)*p(t))/(a*sigma^2);
    D_adj = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*p(t))/(a*sigma^2);
    %kstar_check(t,1) = sum(D<0);
    %kstar_check2(t,1) = sum(D_adj<0);
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
%max(abs(Bind_no - kstar_check))
%max(abs(Bind_no - kstar_check2))
     
%figure(1)
%subplot(1,2,1), hold on, plot(0:20,[p0; p(1:20)],'k'), title('Asset price'), xlabel('Time, t')
%subplot(1,2,2), hold on, plot(1:20,Bind_no(1:20),'k'), title('No. of constrained types'), xlabel('Time, t'), axis([1,inf,-inf,inf])




