function [ combined_p ] = Edgington( p )
% Edgington method of combing p-values when the sum of p-values is less than or equal to 1: Edgington (1972) An additive method for combining probability values from independent experiments
    S = sum(p);
    K = numel(p);
    if(S <= 1)
        combined_p = S^K / factorial(K);
    else
        Fact = factorial(K);
        combined_p = sum(arrayfun(@(j) (-1)^j * nchoosek(K, j) * ( (S-j)^K / Fact ), 0:floor(S)));        
    end    
end

