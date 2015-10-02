function [ selected_features, debug ] = select_features( H, M, residual_Z, annotations, varargin )

    args = inputParser;    
    args.addParamValue('top_perc', 50, @(x) isscalar(x) & x >= 0 & x <= 100); % percent of top ranked cell types we search for maximum seperation point
   
    args.addParamValue('filter_range', false, @(x) islogical(x) ); % Should we limit expression level ?
    args.addParamValue('max_range', 2^16, @(x) isscalar(x) & x >= 0 & x <= 2^16); % max accepted expression
    args.addParamValue('min_range', 0, @(x) isscalar(x) & x >= 0 & x <= 2^16 ); %min accepted expression
    args.addParamValue('pval_threshold', 0.05, @(x) isscalar(x) & x >= 0 & x <= 1); % significane threshold for cuts among top ranked cell-types
    args.addParamValue('fold_threshold', 1, @(x) isscalar(x) & x >= 0 & x <= 1); % Min log2-fold change to consider
    args.addParamValue('epsilon', 1e-1, @(x) isscalar(x) & x >= 0); % epsilon threshold for changes in smoothed condition number
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
        debug.range_exlucded_features = find(G_mask | M_mask);        

%         max_M = max(M, [], 2);
%         max_G = max(G, [], 2);
%         Fold_Z = Modified_zscore(abs(log(max_M ./ max_G)));
%         debug.range_exlucded_features = find(Fold_Z > 1.96);
        
        selected_features = setdiff(selected_features, debug.range_exlucded_features);
        H = H(selected_features, :);    
        G = G(selected_features, :);    
        n = n - numel(debug.range_exlucded_features);
    end
        
%% Pick tissue/celltype selective features
    [~, perm] = sort(G, 2, 'descend');    
    Top_class = num2cell(perm, 2);    
    
    % ********* For debug purposes ***********
    lin_col = perm(:);
    lin_row = repmat((1:n)', q, 1);    
    idx = sub2ind(size(G), lin_row, lin_col);
    G_sorted =  reshape(G(idx), n, q); 
    % ****************************************
    
    fprintf('Computing statistics for top ranked genes ...\n');
    max_col = max(1, floor(size(G, 2)*params.top_perc/100));     
    logPval = ones(n, max_col);
    logFold = zeros(n, max_col);
    for i = 1:max_col     
        fprintf('\tAnalyzing top %d cell-types for each feature ...\n', i);
        Class_members = arrayfun(@(row) H(row, ismember(pure_col_classes, Top_class{row}(1:i))), (1:n)', 'UniformOutput', false);
        Class_stats = cell2mat(cellfun(@(x) [mean(x), var(x)], Class_members, 'UniformOutput', false));
        Class_mean = Class_stats(:, 1); 
        Class_var = Class_stats(:, 2);
        Class_N = cellfun(@(x) numel(x), Class_members);
        Class_var_weighted = Class_var ./ Class_N;
        
        Null_members = arrayfun(@(row) H(row, ismember(pure_col_classes, setdiff(1:q, Top_class{row}(1:i)))), (1:n)', 'UniformOutput', false);
        Null_stats =  cell2mat(cellfun(@(x) [mean(x), var(x)], Null_members, 'UniformOutput', false));
        Null_mean = Null_stats(:, 1); 
        Null_var = Null_stats(:, 2);
        Null_N = l - Class_N;
        Null_var_weighted = Null_var ./ Null_N;
        
        Delta_mean = Class_mean - Null_mean;
        SE = arrayfun(@(x, y) sqrt(x + y), Class_var_weighted, Null_var_weighted);
        T_stat = Delta_mean ./ SE;
        DoF = arrayfun(@(x, n1, y, n2) (x + y)^2 / ( (x^2 / (n1-1)) + (y^2 / (n2-1)) ), Class_var_weighted, Class_N, Null_var_weighted, Null_N); 
        logPval(1:n, i) = -log10(min(1, arrayfun(@(ratio, dfe) tcdf(-ratio, dfe), T_stat, DoF)));
        logFold(1:n, i) = log2(Class_mean ./ Null_mean);
    end

%     cut_size = 5000;
%     [~, perm] = sort(residual_Z, 'descend');
%     max_M = max(M, [], 2);
%     max_G = max(G, [], 2);
    
%     Fold_Z = Modified_zscore(abs(log(max_M ./ max_G)));
%     Max_fold = -Modified_zscore(max(logFold, [], 2));
%     max_pval = -Modified_zscore(max(logPval, [], 2));
%     table = [residual_Z(perm(1:cut_size)), Fold_Z(perm(1:cut_size)), Max_fold(perm(1:cut_size)), max_pval(perm(1:cut_size))];
    

    logFold(logPval < -log10(params.pval_threshold)) = 0;
    feature_score = max(logFold, [], 2);    
    [~, perm] = sort(feature_score, 'descend');
    if(params.fixed_cut == 0)
        potential_cand_count = nnz(params.fold_threshold < feature_score);
        Cond = inf(potential_cand_count-q+1, 1);
        fprintf('Testing condition # over different cuts ...\n');
        for i = q:potential_cand_count 
            G_sub = G(perm(1:i), :);
            Cond(i-q+1) = cond(G_sub);
        end
        Diff_con = abs(diff(Cond));
        Diff_con_smoothed = smooth(Diff_con, 15);
        debug.selectivity_cut_size = find(Diff_con_smoothed < params.epsilon, 1);    
        debug.selectivity_cond = Cond(debug.selectivity_cut_size);
    else
        debug.selectivity_cut_size = nnz(params.fold_threshold < feature_score); % params.fixed_cut;
    end
    debug.selectivity_excluded_features = selected_features(setdiff(1:n, perm(1:debug.selectivity_cut_size)));
    selected_features = selected_features(perm(1:debug.selectivity_cut_size));            
    debug.selected_features = selected_features;          
end

