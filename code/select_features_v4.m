function [ selected_features, debug ] = select_features( H, M, annotations, varargin )

    args = inputParser;    
    args.addParamValue('top_perc', 50, @(x) isscalar(x) & x >= 0 & x <= 100); % percent of top ranked cell types we search for maximum seperation point
   
    args.addParamValue('filter_range', false, @(x) islogical(x) ); % Should we limit expression level ?
    args.addParamValue('max_range', 2^16, @(x) isscalar(x) & x >= 0 & x <= 2^16); % max accepted expression
    args.addParamValue('min_range', 0, @(x) isscalar(x) & x >= 0 & x <= 2^16 ); %min accepted expression

    args.addParamValue('filter_uniformity', false, @(x) islogical(x) ); % Should we filter universally expressed genes first (this is mostly important in reducing computation time as any combination of C will be a solution for uniformly expressed genes)

    args.addParamValue('filter_selectivity', false, @(x) islogical(x) ); % Should we limit expression level ?
    args.addParamValue('pval_threshold', 1e-5, @(x) isscalar(x) & x >= 0 & x <= 1); % significane threshold for cuts among top ranked cell-types
    args.addParamValue('fold_threshold', 2, @(x) isscalar(x) & x >= 0 & x <= 1); % Min log2-fold change to consider
    args.addParamValue('epsilon', 1e-2, @(x) isscalar(x) & x >= 0); % epsilon threshold for changes in smoothed condition number
    args.addParamValue('fixed_cut', 0, @(x) isscalar(x) & x >= 0 & x <= size(H, 1)); % skip cond # and select fixed number of features
    
    args.parse(varargin{:});        
    params = args.Results;   
    
    [n, l] = size(H);
    selected_features = (1:n)';
    [Classes,~ , pure_col_classes] = unique(annotations.pure.Class);
    q = numel(Classes);
    G = cell2mat(arrayfun(@(class) mean(H(:, (pure_col_classes == class)), 2), 1:q, 'UniformOutput', false));
    
%%  Filter out of range probes
    % According to Kawaji et al., small expression values are scewed at
    % lower values, and the effective range of arrays is 3 order of
    % magnitudes. Thus, we remove every feature out of the range [10, 10,000]
    % Ref: Kawaji, H., Lizio, M., Itoh, M., Kanamori-Katayama, M., Kaiho, A., Nishiyori-Sueki, H., … Carninci, P. (2014). Comparison of CAGE and RNA-seq transcriptome profiling using clonally amplified and single-molecule next-generation sequencing. Genome Research, 24(4)
    % Gong et al also suggested the range of 0.5-5K (Gong, T., Hartmann, N., Kohane, I. S., Brinkmann, V., Staedtler, F., Letzkus, M., … Szustakowski, J. D. (2011). Optimal deconvolution of transcriptional profiling data using quadratic programming with application to complex clinical blood samples. PloS One, 6(11), e27156. http://doi.org/10.1371/journal.pone.0027156)
    % Ahn et al. suggested 2^4-2^14 (Ahn, J., Yuan, Y., Parmigiani, G., Suraokar, M. B., Diao, L., Wistuba, I. I., & Wang, W. (2013). DeMix: deconvolution for mixed cancer transcriptomes using raw measured data. Bioinformatics (Oxford, England), 29(15), 1865–71. http://doi.org/10.1093/bioinformatics/btt301)
    if(params.filter_range)
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
        H = H(~exclusion_mask, :);    
        G = G(~exclusion_mask, :);    
        M = M(~exclusion_mask, :); % We don't use M anymore, but for the sake of sticking to main rules of each filter ...
        n = n - nnz(debug.range_exlucded_features); % # of remaining features
    else
       debug.range_exlucded_features = {}; 
    end
%% Filter out uniformly expressed features (Unfinished)
    if(params.filter_uniformity)
        P = normalize(G, 'dim', 2, 'pnorm', 1); % Relative expression of each gene in different tissues
        I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene   
        H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
        [~, nonuniformity_rank] = sort(H);
        Cond_rev = zeros(numel(nonuniformity_rank), 1);
        for i = 0:numel(nonuniformity_rank)-1
            Cond_rev(i) = cond(G(nonuniformity_rank(end-i:end), :));
        end
        Cond = zeros(numel(nonuniformity_rank), 1);
        for i = 1:numel(nonuniformity_rank)
            Cond(i) = cond(G(nonuniformity_rank(1:i), :));
        end        
    end

%% Exclude features that do not have any significant seperation between their highly expressed group and rest of references
    if(params.filter_selectivity)
        Nr = arrayfun(@(c) nnz(pure_col_classes == c), 1:q);
        Perm = zeros(n, q);
        BayesFactors = zeros(n, q);
        for k = 1:q
            fprintf('Selecting markers for %s ..\n', annotations.references.Reference_ID{k});
            tic
            r = Constr_Bayes(H, Nr, k, 1, 1000);
            toc
            r(isinf(r)) = 0;
            [BayesFactors(:, k), Perm(:, k)] = sort(r, 'descend');
        end
        sig_count = sum(BayesFactors > 32);        
        selected_rows =  [];
        Cond = zeros(min(sig_count), 1);
        for i = 1:min(sig_count)
            selected_rows = union(selected_rows, Perm(i, :));
            G_sub = G(selected_rows,:);
            Cond(i) = cond(G_sub);
        end
        
        [~, min_Cond_row] = min(Cond);
        selected_rows = unique(Perm(1:min_Cond_row, :));
        debug.selectivity_excluded = selected_features(setdiff(1:n, selected_rows));
        selected_features = selected_features(selected_rows);
    else
        debug.selectivity_excluded = {};
    end
end

