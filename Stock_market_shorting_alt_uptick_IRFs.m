%Stock market model with short-selling constraint and endogenous shares 
%Simulations with the wealth distribution (Fig. 5)
%Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clear, clc, %close all; 

%Parameter values
H = 1000;
r = 0.1; a = 1;  %benchmark  
betta = 4.5; %3, 4.5
dbar = 0.6; sigma = 1; Zbar = 0.1;  
pf = (dbar - a*sigma^2*Zbar)/r; %Fundamental price
kappa = 0.1;  %Alternative uptick rule

%Coding choices
Iter = 1;  %Iter = 1 turns on iterative algorithm (advisable for large H).
Fixed = 0; %Fixed  = 1: Pick fixed rather than time-varying (fitness-based) population shares. 
n_iter = 6; % no. of iterations (increase for large H)
Back = 0;  %Forward search is the default. Else solve backward from most optimistic type. 
Unconstrained = 0; %Set Unconstrained = 1 to simulate without short-selling constraints. 
T = 30;  %no. of periods

%Preallocate matrices
U = NaN(H,1); Bind = zeros(T,1); Bind_no = Bind;  bound = Bind; D = NaN(H,1); D_lag2 = D; x = NaN(T,1); 
Time = x; Check1 = x; Check11 = x; Demand_vec = NaN(H,T); Wealth_vec = Demand_vec;
count = zeros(T,1); Gini = NaN(T,1); Zero_wealth = Gini; Diff = NaN(H,1); sum_tot = Diff; 

%Dividend shocks 
shock = zeros(T,1);

%Initial values  
p0 = 8; 
x0 = p0 - pf; xlag = p0 - pf; n_init = 1/H*ones(1,H); 
Wealth_init = 50*ones(H,1);

%Disperse beliefs
%b = zeros(H,1); C = b; g = 1.2*ones(H,1); %betta = 3; p0 = 8;
%g(1:H/2) = 0; b(1:H/2) =  linspace(-0.2,0.2,H/2);  C(1:H/2) = 1-abs(b(1:H/2));

%Disperse beliefs 2
rng(10)
b = zeros(H,1); C = b; g = b; g(H/2+1:H) = 1 + 0.4*rand(H/2,1); 
b(1:H/2) =  linspace(-0.2,0.2,H/2); C(1:H/2) = 1-abs(b(1:H/2)); 

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
        %Dem(:,t-2) = NaN(H,1);
        %Compute_fitness_shares_insert
    end
    
        [Beliefs_sort,I] = sort(Beliefs);  
        n_adj = n(I);
       
%Trial unconstrained solution
xstar = n*Beliefs/(1+r);   

if n*Beliefs - min(Beliefs) > a*sigma^2*Zbar && Unconstrained == 0 && cond <=0 
        
        Bind(t) = 1;
        
        %if length(unique(Beliefs)) == H
        %Sort beliefs when there are no ties
        %Check0(t) = sum(n_adj); 
        %else
            %Dum(t) = 1;  
        %    Stock_market_shorting_sort_insert
        %end
            
   %Obtain initial guess for no. short-sellers
   Demand_star = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar)/(a*sigma^2);   
   Demand_star(Demand_star>=0)=0; Demand_star(Demand_star<0)=1; 
   k_init0 = max(sum(Demand_star),1);
   
   %Decide whether to run iterations
   if Iter == 0
       k_init = k_init0;
   else 
       Stock_market_shorting_iterations_insert
       %num_iter(t) = n_iter;
   end  
   
    %Find the equilibrium no.of short-sellers
    %if Back == 0
        for k = k_init:length(Beliefs_sort)-1
            if n_adj(k:end)*Beliefs_sort(k:end) - sum(n_adj(k:end))*Beliefs_sort(k) > a*sigma^2*Zbar && n_adj(k+1:end)*Beliefs_sort(k+1:end) - sum(n_adj(k+1:end))*Beliefs_sort(k+1) <= a*sigma^2*Zbar 
            break
            end
        end
            
    %else 
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

%Check market clearing
    D = (Beliefs + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    D_adj = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*x(t))/(a*sigma^2);
    if Bind(t) == 1
        D(D<0) = 0;  
        D_adj(D_adj<0) = 0; 
    end
    Demand_vec(:,t) = D;
    Check1(t) = abs(n*D - Zbar); 
    Check11(t) = abs(n_adj*D_adj - Zbar);
   
Time(t) = t;

if t==1
Wealth = Wealth_init;
else
Wealth = (1+r)*(Wealth - (x(t-1)+pf)*Demand_vec(:,t-1) ) + ( (x(t)+pf) + dbar + shock(t) )*Demand_vec(:,t-1); 
end

if min(Wealth) < 0
    bound(t) = 1;
end
Wealth(Wealth<0) = realmin;
count(Wealth==realmin) = 1;
Wealth_vec(:,t) = Wealth/max(Wealth);
Wealth_norm = Wealth_vec(:,t);

    %for i=1:length(Wealth)
        %for j=1:length(Wealth)
       
            %Diff(j) = abs(Wealth_norm(i)-Wealth_norm(j));
            
        %end 
        %sum_tot(i) = sum(Diff); 
    %end
    
    %Fast approach to calculate Gini
    V = Wealth_norm';
    Diff_mat = bsxfun(@minus,V(:), V(:).');
    Diff_mat = abs(Diff_mat);
    sum_tot = sum(Diff_mat,1);

    Gini(t) = sum(sum_tot)/(2*(length(Wealth_norm))^2*mean(Wealth_norm));
    
    Zero_wealth(t) = sum(count); 
    Rel_wealth = Wealth/max(Wealth);

end

x = [x0; x]; Time = [0; Time]; 
max(Check1)
max(Check11)
max(bound)

%figure(1)
%hold on, subplot(2,2,4), plot(Time,x,'--k'), hold on, %plot(Bind_no(1:300)/H,'--r')
%subplot(1,2,2), plot(Bind_no(1:200)/max(Bind_no),'--k')
%axis([-inf,inf,-inf,inf]), title('S2: Hetero. fundamentalists, b \in [-0.2,0.2], \beta = 4.5'), ylabel('Price deviation \it{x}'), xlabel('Time')

%figure(2)
%hold on, subplot(2,2,4), plot(Time(2:end),Gini,'--k'), hold on, %plot(Time,Zero_wealth/H,'--r')
%axis([-inf,inf,-inf,inf]), title('S2: Hetero. fundamentalists, b \in [-0.2,0.2], \beta = 4.5'), ylabel('Gini coefficient'), xlabel('Time')
%subplot(1,2,2), plot(Bind_no(1:200)/max(Bind_no),'--k')
%axis([-inf,inf,-inf,inf]),
%subplot(2,2,2), plot(Time,Zero_wealth/H)

figure(3)
subplot(2,3,6), hold on, histogram(Wealth_vec(:,24),40,'FaceColor',[0.5,0.5,0.5]), title('Wealth distribution at t=24'), xlabel('Relative wealth')      
%[0.5,0.5,0.5]

