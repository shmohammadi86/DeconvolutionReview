function [ filtered_rows ] = adaptiveRangeFiltering( M, G, varargin )
    [n, p] = size(M);
    q = size(G, 2);

    args = inputParser;    
    args.addParamValue('filter_type', 'both', @(x) ischar(x)); 
    
    args.parse(varargin{:});        
    params = args.Results;

    X = cell2mat(arrayfun(@(i) log10(sort(smooth(G(:, i)))), 1:q, 'UniformOutput', false));
    
    G_max = max(G, [], 2);
    M_max = max(M, [], 2);
    max_expr = max(G_max, M_max);
    y = log2(sort(max_expr));
    x = 100*(1:size(max_expr, 1)) ./ size(max_expr, 1);        
    mid_point_idx = round(0.5 * size(max_expr, 1));

    min_idx = cut(y(1:mid_point_idx), 'selection_method', 'Elbow', 'selection_half', 'bottom'); min_val_threshold = 2^y(min_idx);    
    max_idx = cut(y(mid_point_idx:end), 'selection_method', 'Elbow', 'selection_half', 'top')+mid_point_idx; max_val_threshold = 2^y(max_idx);  

    max_idx_filter = (max_val_threshold < min(M, [], 2) | max_val_threshold < min(G, [], 2));
    min_idx_filter = (max(M, [], 2) < min_val_threshold | max(G, [], 2) < min_val_threshold);    
    
    switch(params.filter_type)
        case 'both'
            extreme_sample_mask = (1+params.alpha)*Ref_max < y; 
            extreme_ref_mask = (1+params.alpha)*y < Ref_min; 
            bad_feature_mask = extreme_sample_mask | extreme_ref_mask;
        case 'M'
            bad_feature_mask = (1+params.alpha)*Ref_max < y; 
        case 'G'
            bad_feature_mask = (1+params.alpha)*y < Ref_min;
        case 'none'
            bad_feature_mask = false(n, 1);
     end    

end

