%Stock_market_shorting_sort_insert
   
        %Sort beliefs with ties
           [Beliefs_sort,~,rnk] = unique(Beliefs,'sorted');  %Sorting and ranking        
            
            n_adj = NaN(1,max(rnk));
            for h=1:max(1,max(rnk))
                V = zeros(1,H);
                V(rnk == h) = 1;
                n_adj(1,h) = V*transpose(n);
            end
