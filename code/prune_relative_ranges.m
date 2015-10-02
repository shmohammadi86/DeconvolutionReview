function [ bad_feature_mask, G_greaterThan_M_count, M_greaterThan_G_count ] = prune_relative_ranges ( M, G, varargin )
    args = inputParser;    
    args.addParamValue('alpha', 1, @(x) isscalar(x) & x >= 1); 
    args.addParamValue('type', 'median', @(x) ischar(x)); % min/median/max
    
    args.parse(varargin{:});        
    params = args.Results;   
    
    [n, q] = size(G);

    
    P = normalize(G, 'dim', 2, 'pnorm', 1); % Relative expression of each gene in different tissues
    [~, perm] = sort(G, 2, 'descend');             
    lin_col = perm(:);
    lin_row = repmat((1:n)', q, 1);    
    idx = sub2ind(size(G), lin_row, lin_col);

    sorted_refs =  reshape(G(idx), n, q); 
    sorted_probs = reshape(P(idx), n, q); 	
    reverse_sorted_probs = fliplr(sorted_probs);

    min_weighted_avg = params.alpha * sum(reverse_sorted_probs .* sorted_refs, 2);
    max_weighted_avg = (1 / params.alpha) * sum(sorted_probs .* sorted_refs, 2);

    
%     G_greaterThan_M = zeros(n, 1);
%     M_greaterThan_G = zeros(n, 1);
%     for i = 1:q
%         G_greaterThan_M = G_greaterThan_M | (M(:, i) < min_weighted_avg);
%         M_greaterThan_G = M_greaterThan_G  | (max_weighted_avg < M(:, i));
%     end

     [ ~, Median, MAD ] = Modified_zscore(M');
     
%     G_greaterThan_M = Median' + MAD' < min_weighted_avg;
%     M_greaterThan_G = max_weighted_avg < Median' - MAD';

    G_greaterThan_M = Median' < min_weighted_avg;
    M_greaterThan_G = max_weighted_avg < Median';
    
    bad_feature_mask = G_greaterThan_M | M_greaterThan_G;    
    G_greaterThan_M_count = nnz(G_greaterThan_M);
    M_greaterThan_G_count = nnz(M_greaterThan_G);
end
