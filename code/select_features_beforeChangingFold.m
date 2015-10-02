function [ selected_features, debug ] = select_features( H, M, annotations, datasetName, varargin )

    args = inputParser;    
    args.addParamValue('visualize', false, @(x) islogical(x) ); % Should we visualize diagnostic plots?

    args.addParamValue('ntop', 2, @(x) isscalar(x) & x >= 0 & x <= 100); % percent of top ranked cell types we search for maximum seperation point
   
    args.addParamValue('filter_range', false, @(x) islogical(x) ); % Should we limit expression level ?
    args.addParamValue('max_range', 2^16, @(x) isscalar(x) & x >= 0 & x <= 2^16); % max accepted expression
    args.addParamValue('min_range', 0, @(x) isscalar(x) & x >= 0 & x <= 2^16 ); %min accepted expression

    args.addParamValue('filter_selectivity', false, @(x) islogical(x) ); % Should we use ttest to eliminate geatures that are not significantly separated ?
    args.addParamValue('pval_threshold', 1e-2, @(x) isscalar(x) & x >= 0 & x <= 1); % significane threshold for cuts among top ranked cell-types
    
    args.addParamValue('filter_fold', false, @(x) islogical(x) ); % Should we use fold enrichment/ cond # to orthogonalize basis?
    args.addParamValue('fold_threshold', 1, @(x) isscalar(x)); % Min log2-fold change to consider
    args.addParamValue('cond_method', 'Ours', @(x) ischar(x)); % Min log2-fold change to consider

    

    args.addParamValue('filter_uniformity', false, @(x) islogical(x) ); % Should we filter universally expressed genes first (this is mostly important in reducing computation time as any combination of C will be a solution for uniformly expressed genes)
    args.addParamValue('filter_condition', false, @(x) islogical(x) ); % Should we limit expression level ?
    args.addParamValue('epsilon', 1e-2, @(x) isscalar(x) & x >= 0); % epsilon threshold for changes in smoothed condition number
    args.addParamValue('fixed_cut', 0, @(x) isscalar(x) & x >= 0 & x <= size(H, 1)); % skip cond # and select fixed number of features
    
    args.parse(varargin{:});        
    params = args.Results;   
    [n, l] = size(H);
    selected_features = (1:n)';
    [Classes,~ , pure_col_classes] = unique(annotations.pure.Class);
    q = numel(Classes);
    G = cell2mat(arrayfun(@(class) mean(H(:, (pure_col_classes == class)), 2), 1:q, 'UniformOutput', false));
    params.ntop = min(params.ntop, q);
    
    path.plot_output = sprintf('output/sorted/part3/unlimitedRange_%.0e_%d/plots', params.pval_threshold, params.fold_threshold);
    
%%  Filter out of range probes
    % According to Kawaji et al., small expression values are scewed at
    % lower values, and the effective range of arrays is 3 order of
    % magnitudes. Thus, we remove every feature out of the range [10, 10,000]
    % Ref: Kawaji, H., Lizio, M., Itoh, M., Kanamori-Katayama, M., Kaiho, A., Nishiyori-Sueki, H., … Carninci, P. (2014). Comparison of CAGE and RNA-seq transcriptome profiling using clonally amplified and single-molecule next-generation sequencing. Genome Research, 24(4)
    % Gong et al also suggested the range of 0.5-5K (Gong, T., Hartmann, N., Kohane, I. S., Brinkmann, V., Staedtler, F., Letzkus, M., … Szustakowski, J. D. (2011). Optimal deconvolution of transcriptional profiling data using quadratic programming with application to complex clinical blood samples. PloS One, 6(11), e27156. http://doi.org/10.1371/journal.pone.0027156)
    % Ahn et al. suggested 2^4-2^14 (Ahn, J., Yuan, Y., Parmigiani, G., Suraokar, M. B., Diao, L., Wistuba, I. I., & Wang, W. (2013). DeMix: deconvolution for mixed cancer transcriptomes using raw measured data. Bioinformatics (Oxford, England), 29(15), 1865–71. http://doi.org/10.1093/bioinformatics/btt301)
    if(params.filter_range)
        if(numel(selected_features) == 0) 
            return
        end        
        G_mask = min(G, [], 2) < params.min_range | params.max_range < max(G, [], 2);
        M_mask = min(M, [], 2) < params.min_range | params.max_range < max(M, [], 2);
        exclusion_mask = G_mask | M_mask;
        debug.range_exlucded_features = selected_features(exclusion_mask);        

        % Idea of selective high filter
%         max_M = max(M, [], 2);
%         max_G = max(G, [], 2);
%         Fold_Z = Modified_zscore(abs(log(max_M ./ max_G)));
%         debug.range_exlucded_features = find(Fold_Z > 1.96);
        
        selected_features = selected_features(~exclusion_mask);
        debug.selected_features_after_rangeFilter = selected_features;
%         H = H(~exclusion_mask, :);    
%         G = G(~exclusion_mask, :);    
%         M = M(~exclusion_mask, :); % We don't use M anymore, but for the sake of sticking to main rules of each filter ...
        n = n - nnz(debug.range_exlucded_features); % # of remaining features
    else
       debug.range_exlucded_features = {}; 
    end
%% Filter out uniformly expressed features (Unfinished)
    if(params.filter_uniformity)
%         P = normalize(G, 'dim', 2, 'pnorm', 1); % Relative expression of each gene in different tissues
%         I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene   
%         H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
%         [~, nonuniformity_rank] = sort(H);
%         Cond_rev = zeros(numel(nonuniformity_rank), 1);
%         for i = 0:numel(nonuniformity_rank)-1
%             Cond_rev(i) = cond(G(nonuniformity_rank(end-i:end), :));
%         end
%         Cond = zeros(numel(nonuniformity_rank), 1);
%         for i = 1:numel(nonuniformity_rank)
%             Cond(i) = cond(G(nonuniformity_rank(1:i), :));
%         end        
    end
%%
    if(params.filter_selectivity)
        if(numel(selected_features) == 0) 
            return
        end        
        feature_Pvals = ones(n, params.ntop-1);        
        logH = log2(H(selected_features, :));
        ref_mean = cell2mat(arrayfun(@(class) mean(logH(:, (pure_col_classes == class)), 2), 1:q, 'UniformOutput', false));
        ref_var = cell2mat(arrayfun(@(class) var(logH(:, (pure_col_classes == class)), 0, 2), 1:q, 'UniformOutput', false));
        ref_count = repmat(cell2mat(arrayfun(@(class) nnz(pure_col_classes == class), 1:q, 'UniformOutput', false)), n, 1);
        [~, sorted_classes] = sort(ref_mean, 2, 'descend');         
        
        lin_col = sorted_classes(:);
        lin_row = repmat((1:n)', q, 1);    
        idx = sub2ind([n, size(G, 2)], lin_row, lin_col);
        ref_mean_sorted = reshape(ref_mean(idx), n, q);  
        ref_var_sorted  = reshape(ref_var(idx), n, q);  
        ref_count_sorted       = reshape(ref_count(idx), n, q);  

        Class_mean = ref_mean_sorted(:, 1);
        Class_var = ref_var_sorted(:, 1);
        Class_N = ref_count_sorted(:, 1);
        Class_var_weighted = Class_var ./ Class_N;
        for i = 2:params.ntop
            Null_mean = ref_mean_sorted(:, i);
            Null_var = ref_var_sorted(:, i);
            Null_N = ref_count_sorted(:, i);
            Null_var_weighted = Null_var ./ Null_N;   
            
            Delta_mean = Class_mean - Null_mean;
            SE = arrayfun(@(x, y) sqrt(x + y), Class_var_weighted, Null_var_weighted);
            T_stat = Delta_mean ./ SE;
            DoF = arrayfun(@(x, n1, y, n2) (x + y)^2 / ( (x^2 / (n1-1)) + (y^2 / (n2-1)) ), Class_var_weighted, Class_N, Null_var_weighted, Null_N); 
            feature_Pvals(1:n, i-1) = arrayfun(@(ratio, dfe) tcdf(-ratio, dfe), T_stat, DoF);
        end

%         minPval = min(feature_Pvals, [], 2);
%         [~, feature_Qvals] = mafdr(minPval);
%         
%         feature_Qvals = reshape(feature_Qvals, n, params.ntop-1);
%         x = logspace(-1, log10(1/n), 1000);
%         Qval_threshold = x(find(arrayfun(@(tr) (tr * (nnz(min(feature_Qvals, [], 2) <= tr))) < 1, x), 1));
%         significant_features_mask = min(feature_Qvals, [], 2) <= Qval_threshold;

        
%         [~, feature_Qvals] = mafdr(feature_Pvals(:));
%         feature_Qvals = reshape(feature_Qvals, n, params.ntop-1);
%         x = logspace(-1, log10(1/n), 1000);
%         Qval_threshold = x(find(arrayfun(@(tr) (tr * (nnz(min(feature_Qvals, [], 2) <= tr))) < 1, x), 1));
%         significant_features_mask = min(feature_Qvals, [], 2) <= Qval_threshold;

        % Diagnostic plots
        debug.feature_ttset_Pvals = feature_Pvals;
        if(params.visualize)
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);         
            [f, xi] = ksdensity(-log10(feature_Pvals(:)), 'function', 'cdf');
            plot(xi, f, 'LineWidth', 2);
            critical_point = find(0.9 < f, 1);
            f_c = f(critical_point);
            x_c = xi(critical_point);
            x_range = get(gca, 'Xlim');
            y_range = get(gca, 'Ylim');        
            line([x_range(1), x_c], [f_c, f_c], 'LineWidth', 2, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            line([x_c, x_c], [y_range(1), f_c], 'LineWidth', 2, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            set(gca, 'XLim', [0, x_range(2)]);
            set(gca, 'XTick', 1:floor(x_range(2)))                
            xlabel('-log_{10}(p-values)', 'FontSize',12, 'FontWeight','bold');
            ylabel('CDF', 'FontSize',12, 'FontWeight','bold');
            set(gca,'FontSize',12, 'FontWeight','bold');
            export_fig(sprintf('%s/%s_pval_dist.eps', path.plot_output, datasetName), '-eps', '-transparent');    
            close
        end
        best_feature_Pval = min(feature_Pvals, [], 2);
        best_feature_Pvals_corrected = best_feature_Pval;% min(1, best_feature_Pval .* (numel(feature_Pvals)));
        [feature_Pval_corrected_sorted, feature_Pval_corrected_sorted_perm] = sort(best_feature_Pvals_corrected);

        if(params.visualize)
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);        
            G_sub = G(feature_Pval_corrected_sorted_perm(1:1000), :);        
            myclustergram(G_sub, 'ColumnLabels', annotations.references.Reference_ID, 'Pdist', 'correlation', 'Linkage', 'centroid', 'dimension', 1, 'Colormap', colormap(flipud(redgreencmap())), 'OPTIMALLEAFORDER', true, 'ratio', eps);
            title('Top 1K markers after t-test');
            set(gca,'FontSize',12, 'FontWeight','bold');
            export_fig(sprintf('%s/%s_TopPval_clustergram.eps', path.plot_output, datasetName), '-eps', '-transparent');    
            close
        end
        %         [max_val, max_col] = max(G_sub, [], 2);
%         mean_rest = arrayfun(@(r) mean(G_sub(r, setdiff(1:q, max_col(r))), 2), 1:size(G_sub, 1))';
%         max_val = max_val - mean_rest;
%         [~, max_val_perm] = sort(max_val, 'descend');
%         max_col_sorted = max_col(max_val_perm);
%         G_sub = G_sub(max_val_perm, :);        
%         row_perm = cell2mat(arrayfun(@(c) find(max_col_sorted == c), 1:q, 'UniformOutput', false)');
%         imagesc(G_sub(row_perm, :))
%         colormap(flipud(redgreencmap()));
%         set(gca, 'XTick' , 1:q);
%         set(gca, 'XTickLabel', annotations.references.Reference_ID);        
%         set(gca,'FontSize',12, 'FontWeight','bold');


        significant_features_mask = (best_feature_Pvals_corrected  <= max(feature_Pval_corrected_sorted(20*q), params.pval_threshold)); % make sure to at least have 2q features left after selection    
        debug.ttest_excluded = selected_features(~significant_features_mask);
        selected_features = selected_features(significant_features_mask);
        debug.selected_features_after_ttestFilter = selected_features;
        n = n - nnz(debug.ttest_excluded); % # of remaining features
        
    else
        debug.ttest_excluded = {};
    end
%% Filter based on fold ratio
    if(params.filter_fold)
        if(numel(selected_features) == 0) 
            return
        end
            
        logFold = cell2mat(arrayfun(@(c) log2( G(selected_features, c) ./ mean(G(selected_features, setdiff(1:q, c)), 2) ), 1:q, 'UniformOutput', false));
        [sorted_logFold, perm] = sort(logFold, 'descend');        
        
        debug.max_fold = max(logFold, [], 2);
        if(params.visualize);
            Colors = cbrewer('qual', 'Set1', q);            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]); 
            hold all
            xlabel('# features/ cell-type', 'FontSize',12, 'FontWeight','bold');    
            ylabel('log_2(Fold ratio)', 'FontSize',12, 'FontWeight','bold');
            zeroFold_row = find(min(sorted_logFold, [], 2) < 0, 1);
            for col = 1:q
                 plot(1:zeroFold_row, sorted_logFold(1:zeroFold_row, col), 'LineWidth', 2, 'Color', Colors(col, :));        

            end
            set(gca, 'Ylim', [min(min(sorted_logFold(1:zeroFold_row, :))), max(max(sorted_logFold(1:zeroFold_row, :)))]);               
            cut_point = find(min(sorted_logFold, [], 2) < params.fold_threshold, 1);
            line([cut_point, cut_point], get(gca, 'Ylim'), 'LineWidth', 2, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            current_ticks = setdiff(get(gca, 'Xtick'), 0);
            [~, bad_tick] = min(abs(current_ticks - cut_point));
            current_ticks(bad_tick) = [];
            current_ticks = union(current_ticks, cut_point);
            set(gca, 'Xtick', current_ticks);             

            legend(annotations.references.Reference_ID);
            set(gca,'FontSize',12, 'FontWeight','bold');
            export_fig(sprintf('%s/%s_FoldRatio.eps', path.plot_output, datasetName), '-eps', '-transparent');    
            close
        end    

%         idx = 1;
%         Cond = inf(numel(selected_features)*q, 1);
%         current_selected_features = cell(numel(selected_features)*q, 1);
%         marker_diff = inf(numel(selected_features)*q, 1);
%         Gsqr_sum = cell2mat(arrayfun(@(c) cumsum(G(selected_features(perm(:, c)), c).^2), 1:q, 'UniformOutput', false));
%         col_normSqr = zeros(1, q);
%         col_idx = ones(1, q);
%         [~, max_col] = max(sorted_logFold(1, :), [], 2);
%         [~, min_col] = min(sorted_logFold(1, :), [], 2);
%         col_normSqr(max_col) = Gsqr_sum(1, max_col);
%         col_idx(max_col) = 2;
%         current_selected_features{idx} = perm(1, max_col);
%         G_sub = G(selected_features(current_selected_features{idx}), :); % does G_sub'*G_sub look like eye(q) ??
%         Cond(idx) = cond(G_sub);        
%         
%         while(params.fold_threshold < sorted_logFold(col_idx(min_col), min_col))
%            idx = idx + 1;
%            next_cut = find( col_normSqr(max_col)  < Gsqr_sum(:, min_col), 1);
%            col_normSqr(min_col) = Gsqr_sum(next_cut, min_col);
%            
%            current_selected_features{idx} = union(current_selected_features{idx-1}, perm(col_idx(min_col):next_cut, min_col));
%            gene_size_1(idx) = numel(current_selected_features{idx});
%            col_idx(min_col) = next_cut + 1;
%            marker_diff(idx) = mean(diff(sort(col_idx)));
%            
%            [~, min_col] = min(col_normSqr);
%            [~, max_col] = max(col_normSqr);  
%            
%            G_sub = G(selected_features(current_selected_features{idx}), :); % does G_sub'*G_sub look like eye(q) ??
%            Cond(idx) = cond(G_sub);
%         end
%         Cond(idx+1:end) = []; Cond(1) = [];
%         current_selected_features(idx+1:end) = []; current_selected_features(1) = []; 
%         marker_diff(idx+1:end) = []; marker_diff(1:2) = [];
%         
%         [~, min_cond_idx] = min(smooth(Cond, 15));
        
        max_row = min(numel(selected_features), max(200, find(min(sorted_logFold, [], 2) < params.fold_threshold, 1)));
        min_row = 10;% min(50, numel(selected_features));
        Cond2 = inf(max_row, 1);
        for i = min_row:max_row
            current_features = unique(perm(1:i, :));
            gene_size_2(i) = numel(current_features);
            G_sub= G(selected_features(current_features, :), :);
            Cond2(i-min_row+1) = cond(G_sub);
        end
        
        smoothed_cond = smooth(Cond2, 13);
        if(strcmp(params.cond_method, 'Newman'))
            [~, min_idx] = min(smoothed_cond);
            best_features = selected_features(unique(perm(1:(min_idx+min_row), :)));
        elseif(strcmp(params.cond_method, 'Ours'))
            [~, min_idx] = min(smooth(Cond, 13));
            best_features = selected_features(current_selected_features{min_idx});            
        else
            error('Unknown method %s to find min cond #', params.cond_method);
        end

        if(params.visualize)
%             % Comparison of condition selection methods   
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
            x_ticks = q.*(min_row - 1 + (1:min(numel(Cond2), 2*min_idx)));
            plot(x_ticks,smoothed_cond(1:min(numel(Cond2), 2*min_idx)), 'LineWidth', 2, 'Color', Colors(1, :));
            xlabel('# features', 'FontSize',12, 'FontWeight','bold');    
            ylabel('Condition #', 'FontSize',12, 'FontWeight','bold');
            set(gca,'FontSize',12, 'FontWeight','bold');        
            
%             semilogx(gene_size_2, smooth(Cond2, 13), 'LineWidth', 2, 'Color', Colors(1, :))
%             hold all
%             loglog(gene_size_1(2:end), smooth(Cond, 13),'LineWidth', 2, 'Color', Colors(2, :));        
%             xlabel('log_{10} (# features)', 'FontSize',12, 'FontWeight','bold');    
%             ylabel('log_{10}(Condition #)', 'FontSize',12, 'FontWeight','bold');
%             legend('Newman', 'Ours');
%             grid on
            export_fig(sprintf('%s/%s_CondNum.eps', path.plot_output, datasetName), '-eps', '-transparent');    
            close
        end;
        
        debug.fold_excluded = setdiff(selected_features, best_features);
        selected_features = best_features;
        debug.selected_features_after_foldFilter = selected_features;  
        n = n - nnz(debug.fold_excluded); % # of remaining features
    else
        debug.fold_excluded = {};        
    end
    
%% Visualize final basis
    if(params.visualize)
        set(gcf,'units','normalized','outerposition',[0 0 1 1]); 
        myclustergram(G(selected_features, :), 'ColumnLabels', annotations.references.Reference_ID, 'Pdist', 'cosine', 'Linkage', 'average', 'dimension', 1, 'Colormap', colormap(flipud(redgreencmap())), 'OPTIMALLEAFORDER', true, 'ratio', eps)
        set(gca, 'XTick', 1:size(G, 2));
        set(gca, 'XTickLabel', annotations.references.Reference_ID);
        set(gca,'FontSize',12, 'FontWeight','bold')
        title(sprintf('Total number of marker genes = %d', numel(selected_features)),'FontSize',14);
        export_fig(sprintf('%s/%s_basis.eps', path.plot_output, datasetName), '-eps', '-transparent', '-nocrop');
        close
    end
end

