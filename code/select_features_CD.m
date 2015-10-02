function [ selected_features, debug ] = select_features( H, M, annotations, datasetName, varargin )

    args = inputParser;    
    args.addParamValue('visualize', false, @(x) islogical(x) ); % Should we visualize diagnostic plots?
    args.addParamValue('plot_path', '', @(x) ischar(x) ); % in case of visualization, where should we store the plots?

    args.addParamValue('ntop', 2, @(x) isscalar(x) & x >= 0 & x <= 100); % percent of top ranked cell types we search for maximum seperation point
   
    args.addParamValue('filter_range', false, @(x) islogical(x) ); % Should we limit expression level ?
    args.addParamValue('max_range', 2^16, @(x) isscalar(x) & x >= 0 & x <= 2^16); % max accepted expression
    args.addParamValue('min_range', 0, @(x) isscalar(x) & x >= 0 & x <= 2^16 ); %min accepted expression

    args.addParamValue('filter_selectivity', false, @(x) islogical(x) ); % Should we use ttest to eliminate geatures that are not significantly separated ?
    args.addParamValue('pval_threshold', 1e-3, @(x) isscalar(x) & x >= 0 & x <= 1); % significane threshold for cuts among top ranked cell-types
    
    args.addParamValue('filter_fold', false, @(x) islogical(x) ); % Should we use fold enrichment/ cond # to orthogonalize basis?
    args.addParamValue('fold_threshold', 2, @(x) isscalar(x)); % Min log2-fold change to consider
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
    
    if( params.visualize == true && isempty(params.plot_path) )
        warning('select_features:plot_path', 'Visualization is disabled due to lack of "plot_path" parameter');
        params.visualize = false;
    end

    debug.range_exlucded_features = {};
    debug.ttest_excluded = {};
    debug.fold_excluded = {};
     
       
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
        
        selected_features = selected_features(~exclusion_mask);
        debug.selected_features_after_rangeFilter = selected_features;
        n = n - nnz(debug.range_exlucded_features); % # of remaining features
    else
       debug.range_exlucded_features = {}; 
    end
%% Filter selectivity based on t-test    
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

        min_pval = min(feature_Pvals, [], 2);        
%         [~, min_qval] = mafdr(min_pval);
%         min_qval = min(1, min_pval* numel(selected_features));
        min_qval = min_pval;

        debug.feature_ttset_Pvals = min_pval;
        debug.feature_ttset_Qvals = min_qval;
              
        % Finde the threshold over min_qval such that the expected number
        % of FP among Positives is less than 1
%         x = logspace(log10(max(min_qval)), log10(min(min_qval)), 1000);
%         Qval_threshold = x(find(arrayfun(@(tr) (tr * (nnz(min_qval <= tr))) < 1, x), 1));
        significant_features_mask = min_qval <= params.pval_threshold; % Qval_threshold;
   
        
        % Plot after correction
        if(params.visualize)
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);                 
%             loglog((sort(min_pval)), 'r', 'LineWidth', 2);
%             hold all
%             loglog((sort(min_qval)), 'g', 'LineWidth', 2);
%             set(gca,'FontSize',10, 'FontWeight','bold');           
%             
%             xlabel('log_{10}(# features)','FontSize',12, 'FontWeight','bold');
%             ylabel('log_{10}(p-value)','FontSize',12, 'FontWeight','bold');
%             legend({'p-value', 'corrected p-value'},'FontSize',14, 'FontWeight','bold');
%             export_fig(sprintf('%s/%s_pval_plot.eps', params.plot_path, datasetName), '-eps', '-transparent');    
%             close            
%             
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);     
%             pts = linspace(0, max(max(-log10(min_pval)), max(-log10(min_qval))), 200);
%             [f, xi] = ksdensity(-log10(min_pval), pts, 'function', 'pdf');
%             [f2, xi2] = ksdensity(-log10(min_qval), pts, 'function', 'pdf');
%             plot(xi, f, xi2, f2, 'LineWidth', 2);
%             max_val = -log10(min(min(min_pval, min_qval)));
%             set(gca, 'XLim', [0, max_val]);
%             set(gca, 'XTick', 1:floor(max_val))                
%             line([-log10(params.pval_threshold),-log10(params.pval_threshold)], get(gca, 'YLim'), 'LineWidth', 2, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
%             xlabel('-log_{10}(p/q-values)', 'FontSize',12, 'FontWeight','bold');
%             ylabel('Density', 'FontSize',12, 'FontWeight','bold');
%             legend({'p-value', 'q-value'}, 'FontSize',14, 'FontWeight','bold');
%             set(gca,'FontSize',12, 'FontWeight','bold');           
%             export_fig(sprintf('%s/%s_pval_dist.eps', params.plot_path, datasetName), '-eps', '-transparent');    
%             close
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);     
            pts = linspace(0, max(max(-log10(min_pval)), max(-log10(min_qval))), 200);
            [f, xi] = ksdensity(-log10(min_pval), pts, 'function', 'pdf');
            plot(xi, f, 'LineWidth', 2);
            max_val = -log10(min(min(min_pval, min_qval)));
            set(gca, 'XLim', [0, max_val]);
            set(gca, 'XTick', 1:floor(max_val))                
            line([-log10(params.pval_threshold),-log10(params.pval_threshold)], get(gca, 'YLim'), 'LineWidth', 2, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            xlabel('-log_{10}(p-values)', 'FontSize',12, 'FontWeight','bold');
            ylabel('Density', 'FontSize',12, 'FontWeight','bold');
            legend({'p-value'}, 'FontSize',14, 'FontWeight','bold');
            set(gca,'FontSize',12, 'FontWeight','bold');           
            export_fig(sprintf('%s/%s_pval_dist.eps', params.plot_path, datasetName), '-eps', '-transparent');    
            close            
        end
        
        [~, qval_perm] = sort(min_qval);
        qval_sorted_features = selected_features(qval_perm);
        G_sub = G(qval_sorted_features(1:500), :);
        G_sub_ratio = cell2mat(arrayfun(@(c) G_sub(:, c) ./ mean(G_sub(:, setdiff(1:q, c)), 2) , 1:q, 'UniformOutput', false));
        if(params.visualize)
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);                            
            myclustergram(G_sub_ratio, 'ColumnLabels', annotations.references.Reference_ID, 'Colormap', colormap(flipud(redgreencmap())), 'ratio', eps, 'OPTIMALLEAFORDER', true, 'Pdist', 'correlation');%, 'Pdist', 'correlation', 'OPTIMALLEAFORDER', true, 'ratio', eps);            
            set(gca,'FontSize',12, 'FontWeight','bold');
            title('Top 500 markers after t-test', 'FontSize',14, 'FontWeight','bold');
            export_fig(sprintf('%s/%s_Top500Qval_clustergram.eps', params.plot_path, datasetName), '-eps', '-transparent');    
            close
        end

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
        
        G_sub = G(selected_features, :);
        [G_sub_sorted, mean_perm] = sort(G_sub, 2, 'descend');       
        
        logFold = zeros(numel(selected_features), params.ntop-1);
        for i = 2:params.ntop
            logFold(:, i-1) = log2(G_sub_sorted(:, 1) ./ G_sub_sorted(:, i));
        end        
        debug.max_fold = max(logFold, [], 2);
        
        if(params.visualize);
            Colors = cbrewer('qual', 'Set1', q);            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]); 
            hold all
            
            x_range = [min(debug.max_fold), max(debug.max_fold)];
            pts = linspace(x_range(1), x_range(2), 200);
            for col = 1:q
                [f, xi] = ksdensity(debug.max_fold(mean_perm(:, 1) == col), pts, 'function', 'pdf');
                plot(xi, f, 'LineWidth', 2, 'Color', Colors(col, :));
            end
                        
            set(gca,'FontSize',10, 'FontWeight','bold');
            xlabel('log_2(FoldRatio) for marker genes', 'FontSize',12, 'FontWeight','bold');    
            ylabel('Density', 'FontSize',12, 'FontWeight','bold');
            line([params.fold_threshold, params.fold_threshold], get(gca, 'ylim'), 'LineWidth', 2, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            legend(annotations.references.Reference_ID, 'FontSize',14, 'FontWeight','bold');            
            export_fig(sprintf('%s/%s_FoldRatio.eps', params.plot_path, datasetName), '-eps', '-transparent');    
            close
        end    


        [sorted_logFold, sorted_logFold_perm] = sort(debug.max_fold, 'descend');
        max_row = min(find(min(sorted_logFold, [], 2) < params.fold_threshold, 1) - 1);
%         min_row = 10*q;
%         Cond = inf(max_row - min_row + 1, 1);
%         if(numel(Cond) == 0)
%             warning('select_feature:bad_cond', 'No feature left to assess condition #');
%             return
%         end
%         for i = min_row:max_row
%             current_features = sorted_logFold_perm(1:i, :);
%             G_sub= G(selected_features(current_features, :), :);
%             Cond(i-min_row+1) = cond(G_sub);
%         end
%         smoothed_cond = smooth(Cond, 13);
%         [~, min_idx] = min(smoothed_cond);
        
%         if(params.visualize)            
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
%             plot_range = 1:min(numel(Cond), 2*min_idx);
%             plot(q:plot_range(end)+q-1,smoothed_cond(plot_range), 'LineWidth', 2, 'Color', Colors(1, :));
%             xlabel('# features', 'FontSize',12, 'FontWeight','bold');    
%             ylabel('Condition #', 'FontSize',12, 'FontWeight','bold');
%             set(gca,'FontSize',12, 'FontWeight','bold');        
%             set(gca, 'xtick', 0:10*q:numel(Cond));
%             set(gca, 'xlim', [q, plot_range(end)+q]);
%             export_fig(sprintf('%s/%s_CondNum.eps', params.plot_path, datasetName), '-eps', '-transparent');    
%             close
%         end;
%         
%         best_features = selected_features(sorted_logFold_perm(1:(min_idx+min_row)), :);
        best_features = selected_features(sorted_logFold_perm(1:max_row));
        debug.fold_excluded = setdiff(selected_features, best_features);
        selected_features = best_features;
        debug.selected_features_after_foldFilter = selected_features;  
        n = n - nnz(debug.fold_excluded); % # of remaining features
    else
        debug.fold_excluded = {};        
    end
    
%% Visualize final basis
%     if(params.visualize)
%         G_sub = G(selected_features, :);
%         G_sub_ratio = cell2mat(arrayfun(@(c) G_sub(:, c) ./ mean(G_sub(:, setdiff(1:q, c)), 2) , 1:q, 'UniformOutput', false));        
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]); 
%         myclustergram(G_sub_ratio, 'ColumnLabels', annotations.references.Reference_ID, 'Pdist', 'cosine', 'Linkage', 'average', 'dimension', 1, 'Colormap', colormap(flipud(redgreencmap())), 'OPTIMALLEAFORDER', true, 'ratio', eps)
%         set(gca, 'XTick', 1:size(G, 2));
%         set(gca, 'XTickLabel', annotations.references.Reference_ID);
%         set(gca,'FontSize',12, 'FontWeight','bold')
%         title(sprintf('Total number of marker genes = %d', numel(selected_features)),'FontSize',14);
%         export_fig(sprintf('%s/%s_basis.eps', params.plot_path, datasetName), '-eps', '-transparent', '-nocrop');
%         close
%     end
end

