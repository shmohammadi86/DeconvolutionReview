function [overall_score, sample_pvals] = evaluateC_pval(known_C, estimated_C)

    % JUST to make sure everything is absolutely OKAY
    known_C(known_C < 0) = 0;
    known_C = normalize(known_C, 'pnorm', 1, 'dim', 1); 
    estimated_C(estimated_C < 0) = 0;
    estimated_C = normalize(estimated_C, 'pnorm', 1, 'dim', 1); 
    
    sample_pvals = ones(size(known_C, 2), 1); % for each sample
    for i = 1:size(known_C, 2)
        [~, sample_pvals(i)] = kstest2(known_C(:, i), estimated_C(:, i)); %null hypothesis that the data in vectors x1 and x2 are from the same continuous distribution        
    end
    sample_pvals = 1-sample_pvals; % p-value of difference
    overall_score = -log10(Edgington(sample_pvals));        
end
    
