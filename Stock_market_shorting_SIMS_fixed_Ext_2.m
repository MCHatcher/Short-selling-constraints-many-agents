%Stock market model with short-selling constraint and endogenous shares: simulations 
%Last updated: Dec 12, 2022. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clear, clc, %close all; 

%------------------
%Parameter values
%------------------
H = 3000;
r = 0.1; a = 1; 
dbar = 1.1; Zbar = 0.1;  
cbar = 1;
r = r-cbar; %Adjustment for demand function

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
%Set initial random seed
rng(1)
var_h = a*( 0.9 + 0.2*rand(1,H) );
rng(5)  %Set random seed
a_tild = 1./var_h;  a_tild_orig = a_tild;
a_tild_prime = a_tild'; a_tild_orig_prime = a_tild_orig'; 

Bind = zeros(T,1); Bind_no = Bind;  D = NaN(H,1); D_lag2 = D; p = NaN(T,1); Time = p;
Check1 = NaN(T,1); Check11 = Check1; Beliefs = NaN(H,1); Demand_star = D; Demand_star1 = D;

%--------------------------
%Generate dividend shocks 
%--------------------------
%rng(1), sigma_d  = 0.0099;   
%pd = makedist('Normal','mu',0,'sigma',sigma_d);  %Truncated normal distribution
%pd_t = truncate(pd,-dbar,dbar);
%rng(1), shock = random(pd_t,T,1);   
%NB. For timed sims, store shock in workspace and comment out other lines

%-------------------------------
%Initial values and predictors 
%-------------------------------
n_init = 1/H*ones(1,H); 
n = n_init;  n_tild = a_tild.*n;
Zbar_h = Zbar*1./a_tild;  
Zbar_h_orig = Zbar_h'; 

%Fundamental price
pf = (dbar - Zbar/sum(n_tild))/(r+cbar); %Fundamental price
p0 = pf + 0.6; plag = p0; plag1 = p0; plag2 = p0;

%Disperse beliefs
b = zeros(H,1);  g1 = 0.5*rand(H,1); g2 = 0.2*rand(H,1); 
g3 = -0.1*rand(H,1); g4 = -0.1*rand(H,1);
gf = -0.2 -0.6*rand(H,1); u_stack = 0.04*randn(H,T); u_stack(1:H,1:10) = 0;

%Trend followers
g3(1:H/3) = 0; g4(1:H/3) = 0; gf(1:H/3) = 0; %u_stack(1:H/3,1:T) = 0;

%Contrarians
g1(H/3+1:H*2/3) = 0; g2(H/3+1:H*2/3) = 0; gf(H/3+1:H*2/3) = 0; %u_stack(H/3+1:H*2/3,1:T) = 0;

%Fundamentalists
g1(H*2/3+1:end) = 0; g2(H*2/3+1:end) = 0; g3(H*2/3+1:end) = 0; g4(H*2/3+1:end) = 0;


for t=1:T 
    
    %Beliefs = NaN(H,1);
    
    u = u_stack(:,t);
    
    if t==1
        Beliefs = b + (g1+g3)*(p0 - plag) + (g2+g4)*(plag-plag2) + gf*(p0-pf) + u + dbar -Zbar_h_orig; 
    elseif t==2
        Beliefs = b + (g1+g3)*(p(t-1) - p0) + (g2+g4)*(p0 -plag) + gf*(p(t-1)-pf) + u + dbar -Zbar_h_orig;
    elseif t==3
        Beliefs = b + (g1+g3)*(p(t-1) - p(t-2)) + (g2+g4)*(p(t-2) -p0) + gf*(p(t-1)-pf) + u + dbar -Zbar_h_orig;
    else
        Beliefs = b + (g1+g3)*(p(t-1) - p(t-2)) + (g2+g4)*(p(t-2) -p(t-3)) + gf*(p(t-1)-pf) + u + dbar -Zbar_h_orig;
    end
 
       
%------------------------------    
%Trial unconstrained solution
%------------------------------
pstar = n_tild*Beliefs/( (1+r)*sum(n_tild) ); 

Beliefs_adj = Beliefs + Zbar_h_orig;

        [Beliefs_sort,I] = sort(Beliefs_adj);  
        n_adj = n(I);  n_adj_tild = n_tild(I);  a_tild_adj =  a_tild(I);
        a_tild_prime = a_tild_adj';
        Beliefs_sort = Beliefs(I);  Beliefs_adj = Beliefs_adj(I);
        Zbar_h = Zbar_h_orig(I);  Zbar_h_prime = Zbar_h';

if n_tild*Beliefs - sum(n_tild)*min(Beliefs_adj) > 0 && Unconstrained == 0 
        
        Bind(t) = 1;
        
%Sort beliefs when there are ties (uncomment to use, not essential)
        %if length(unique(Beliefs)) ~= H
        %    run Stock_market_shorting_sort_insert
        %end
    
%--------------------------------------------       
%Obtain initial guess for no. short-sellers
%--------------------------------------------
   Demand_star = a_tild_prime.*(Beliefs_sort + Zbar_h - (1+r)*pstar);    
   k_init0 = sum(Demand_star<0);
   
   %run Stock_market_shorting_iterations_insert
   Stock_market_shorting_iterations_insert_Ext_2 
   
%-----------------------------------------   
%Find the equilibrium no.of short-sellers
%-----------------------------------------
        for k = k_init:length(Beliefs_sort)-1
            if n_adj_tild(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj_tild(k+1:end))*Beliefs_adj(k) > (1-sum(n_adj(k+1:end)))*Zbar && n_adj_tild(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj_tild(k+1:end))*Beliefs_adj(k+1) <= (1-sum(n_adj(k+1:end)))*Zbar 
            break
            end
        end

kstar = k;  Bind_no(t) = k;   %No. of constrained types
          
       p(t) = ( n_adj_tild(kstar+1:end)*Beliefs_sort(kstar+1:end) - sum(n_adj(1:kstar))*Zbar  ) / ( (1+r)*sum(n_adj_tild(kstar+1:end)) );   
       
else 
        p(t) = pstar;   %Solution when SS constraints are slack           
end

%-------------------------------------------------------
%Check market clearing (uncomment for acccuracy checks)
%-------------------------------------------------------
    D = a_tild_orig_prime.*(Beliefs + Zbar_h_orig - (1+r)*p(t));
    D_adj = a_tild_prime.*(Beliefs_sort + Zbar_h - (1+r)*p(t));
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
     
figure(1)
subplot(2,2,3), hold on, plot(0:20,[p0; p(1:20)],'k'), title('Asset price'), xlabel('Time, t')
subplot(2,2,4), hold on, plot(1:20,Bind_no(1:20),'k'), title('No. of constrained types'), xlabel('Time, t')




