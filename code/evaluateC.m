function [KL, mAD, RMSD, R2D] = evaluateC(known_C, estimated_C, sample_no)
    function [KL, mAD, RMSD, R2D]  = computeStats(C, C_hat)
        C_hat(isnan(C_hat)) = 0;
        perc = 100*C(:); perc_hat = 100*C_hat(:);
        
        p = perc ./ sum(perc);  p(p < eps) = eps; p = p ./ sum(p);        
        q = perc_hat ./ sum(perc_hat); q(q < eps) = eps; q = q ./ sum(q);
        KL = sum(p.*log(p./q));

        absErr = abs(perc - perc_hat);
        mAD = mean(absErr);
        RMSD = sqrt(mean(absErr .^ 2));
%         sample_pval = ones(size(C, 2), 1);
%         for i = 1:size(C, 2)
%             [~, sample_pval(i)] = kstest2(C(:, i), C_hat(:, i)); %null hypothesis that the data in vectors x1 and x2 are from the same continuous distribution        
%         end
%         R2D = -log10(Edgington(1-sample_pvals));        
        
%         [~, pvals] = corr([C, C_hat]); % problem when there are only 2 celltypes (zero variance)
%         sample_pvals = diag(pvals(1:size(C, 2), size(C, 2)+1:2*size(C, 2)));
%         R2D = -log10(Edgington(sample_pvals));        
        R2D = 1 - corr(perc, perc_hat);
    end

    % JUST to make sure everything is absolutely OKAY
    known_C(known_C < 0) = 0;
    known_C = normalize(known_C, 'pnorm', 1, 'dim', 1); 
    estimated_C(estimated_C < 0) = 0;
    estimated_C = normalize(estimated_C, 'pnorm', 1, 'dim', 1);     
    
%     rand_C = reshape(cell2mat(arrayfun(@(i) normalize(rand(size(known_C)), 'pnorm', 1, 'dim', 1), 1:sample_no, 'UniformOutput', false)), size(known_C, 1), size(known_C, 2), sample_no);
%     rand_sample_mAD = cell2mat(arrayfun(@(i) mean(abs(known_C - rand_C(:, :, i)), 1), 1:sample_no, 'UniformOutput', false)');
%     
%     for i = 1:sample_no
%         rand_C = normalize(rand(size(known_C)), 'pnorm', 1, 'dim', 1)
%     end
    [KL, mAD, RMSD, R2D] = computeStats(known_C, estimated_C);

end
    
