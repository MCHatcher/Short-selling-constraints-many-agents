%Stock_market_shorting_sort_insert
%This code allows the sorting of beliefs and pop. shares in case of ties
%The main code works without this part
%Last updated: July 2022. Written by Michael Hatcher (m.c.hatcher@soton.ac.uk) 
   
        %Sort beliefs with ties
           [Beliefs_sort,~,rnk] = unique(Beliefs,'sorted');  %Sorting and ranking        
            
            n_adj = NaN(1,max(rnk));
            for h=1:max(1,max(rnk))
                V_00 = zeros(1,H);
                V_00(rnk == h) = 1;
                n_adj(1,h) = V_00*transpose(n);
            end
