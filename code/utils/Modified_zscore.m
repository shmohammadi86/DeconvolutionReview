function [ Z, Median, MAD ] = Modified_zscore( x )
    if(isrow(x))
        x = x';
    end
    Median = median(x);
    Delta = bsxfun(@minus, x,  Median);    
    MAD = median( abs(Delta));

    Z = bsxfun(@rdivide, 0.6745*(Delta), MAD);
end

