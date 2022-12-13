%Stock market iterations (Dec 9, 2022)
%Iterative prcoedure to be used in conjunction with Stock_market_shorting.m. The file is run
%automatically if Iter = 1. It is advisable to set Iter = 1 when H (no. of types) is large. 

j = k_init0;

for i=1:n_iter
   
   xstar1 = ( n_adj(j+1:end)*Beliefs_sort(j+1:end) - sum(n_adj(1:j))*a*sigma^2*Zbar  ) / ( sum(n_adj_tild(j+1:end)) );     
   
   Demand_star1 = (Beliefs_sort + a*sigma^2*Zbar - R_sort*xstar1)/(a*sigma^2);
   
 if j == sum(Demand_star1<0)
     break
 end
  
  j = max(sum(Demand_star1<0),1);  
   
end

k_init = j;
   
 