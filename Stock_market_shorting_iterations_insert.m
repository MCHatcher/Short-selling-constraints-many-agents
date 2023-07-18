%Use of iterations to speed up the solution with short-selling constraints
%The file is run automatically if Iter = 1. It is advisable to set Iter = 1 when H (no. of types) is large.
%Last updated: July 26, 2022. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

j = k_init0;

for i=1:n_iter
   
   xstar1 = ( n_adj(j+1:end)*Beliefs_sort(j+1:end) - sum(n_adj(1:j))*a*sigma^2*Zbar  ) / ( (1+r)*sum(n_adj(j+1:end)) );     
   
   Demand_star1 = (Beliefs_sort + a*sigma^2*Zbar - (1+r)*xstar1)/(a*sigma^2);
   
 if j == sum(Demand_star1<0)
     break
 end
  
  j = max(sum(Demand_star1<0),1);  
   
end

k_init = j;
   
 