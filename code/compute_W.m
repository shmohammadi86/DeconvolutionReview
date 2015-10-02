function [ cc ] = compute_W( H, G, M, knownC, annotations )
    [Classes,~ , pure_col_classes] = unique(annotations.pure.Class);
    q = numel(Classes);            
    ref_col_range = 1:q;
 
    out_of_range_features = find(min(H, [], 2) < 2^4 | 2^16 < max(H, [], 2));
    selected_features = setdiff(1:size(H, 1), out_of_range_features);
    %H = H(selected_features, :);

    [n_sub, l] = size(H);    

    V = zeros(n_sub, q);
    for class = ref_col_range
        class_mask = (pure_col_classes == class);
        V(:, class) = var(H, 0, 2);            
    end
    
    
%     for i = 1:size(M, 2)
%     	w = sum(V, 2) ./ (q^2);
%  		cc(:, i) = lscov(M(:, i), G)%, w); 
%     end
 	w = sum(V, 2) ./ (q^2);
%     cc =(bsxfun(@times, w, bsxfun(@ldivide, sum(M.^2), M))'*G)';

    cc = (bsxfun(@ldivide, sum(M.^2), M)'*bsxfun(@ldivide, sum(G.^2), G))';
    cc = (M'*G)';
    cc(cc < 0) = 0;
	cc = bsxfun(@ldivide, sum(cc), cc); 

end

