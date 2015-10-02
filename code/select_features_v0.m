function [ final_features, debug ] = select_features( H, annotations, varargin )
%% Can be called multiple times to prune iteratively based on different criteria

    debug.opt = -1; % For Huber/SVR, it holds the optimal hyper-parameter for the loss function

    [Classes,~ , pure_col_classes] = unique(annotations.pure.Class);
    q = numel(Classes);            
    ref_col_range = 1:q;

    % According to Kawaji et al., small expression values are scewed at
    % lower values, and the effective range of arrays is 3 order of
    % magnitudes. Thus, we remove every feature out of the range [10, 10,000]
    % Ref: Kawaji, H., Lizio, M., Itoh, M., Kanamori-Katayama, M., Kaiho, A., Nishiyori-Sueki, H., … Carninci, P. (2014). Comparison of CAGE and RNA-seq transcriptome profiling using clonally amplified and single-molecule next-generation sequencing. Genome Research, 24(4)
    out_of_range_features = find(min(H, [], 2) < 10 | 10000 < max(H, [], 2));
    selected_features = setdiff(1:size(H, 1), out_of_range_features);
    H = H(selected_features, :);

    [n_sub, l] = size(H);    
    selected_features_range = 1:n_sub;    
    pure_col_range = 1:l;
    
    
    args = inputParser;    
    args.addParamValue('ranking_method', 'Diff', @(x) ischar(x)); % Options) Diff: differential expression, Var: Biological variation in mixtures, Uni: Uniformity of expression (HK/ unexpressed)
    args.addParamValue('partition_method', 'max_gap', @(x) ischar(x)); % Options) group_vs_rest: compares each group against all other groups, max_gap: Pairwise comparison of each group versus other groups, and select the maximum gap in diff expression of gene in group vs any other group
    args.addParamValue('selection_method', 'eq', @(x) ischar(x)); % Options) eq: selects the same number of genes/ cell-type such that # features = q*G, where q is # celltypes and G is # features/celltype , neq: selects # feature from all genes in all cell types

    args.addParamValue('min_features', q, @(x) isscalar(x) & x >= 0 & x <= size(H, 1)); % We need at least ~one marker/ celltype
    args.addParamValue('max_features', 500*q, @(x) isscalar(x) & x >= 0 & x <= size(H, 1)); % ~ 500 markers/ celltype
    args.parse(varargin{:});        
    params = args.Results;    
        

    params.max_features = min(params.max_features, n_sub); % We can't select more features than the total we have!
    
    G = zeros(n_sub, q);
    switch params.ranking_method
        case 'Diff'
            Stats = zeros(n_sub, q);           
            
            for class = ref_col_range
                
                class_mask = (pure_col_classes == class);
                G(:, class) = mean(H(:, class_mask), 2);
                
                other_classes = setdiff(ref_col_range, class);
                switch params.partition_method
                    case 'group_vs_rest'
                        ctrl_mask = ismember(pure_col_classes, other_classes);
                        curr_stats = -inf(n_sub, 1);
                        case_mat = log2(H(:, class_mask));
                        ctrl_mat = log2(H(:, ctrl_mask));
                        noVar_features = (var(case_mat, 0, 2) < 1e-5) | (var(ctrl_mat, 0, 2) < 1e-5);
                        stats = chdir(ctrl_mat(~noVar_features, :), case_mat(~noVar_features, :), 1);
                        curr_stats(~noVar_features) = stats;
                        
                    case 'max_gap'
                        best_stats = -inf(n_sub , 1);
                        for other_class = other_classes
                            ctrl_mask = (pure_col_classes == other_class);

                            curr_stats = -inf(n_sub, 1);
                            case_mat = log2(H(:, class_mask));
                            ctrl_mat = log2(H(:, ctrl_mask));
                            noVar_features = (var(case_mat, 0, 2) < 1e-5) | (var(ctrl_mat, 0, 2) < 1e-5);
                            good_features = selected_features(~noVar_features);
                            stats = chdir(case_mat(good_features, :), case_mat(ctrl_features, :), 1);
                            curr_stats(~noVar_features) = stats;
                                                    
                            best_stats = max(best_stats, cur_stats);
                        end
                        cur_stats = best_stats;
                        
                    otherwise
                        error('Unknown partition method %s\n_sub', params.papartition_method);
                end
                Stats(:, class) = curr_stats;
            end
            switch params.selection_method
                case 'eq'
                    [~, Stats_perm] = sort(Stats, 'descend');
                    sorted_rows = selected_features_range(Stats_perm);  
                    
                    % check to makes sure all cell types have at least
                    % params.min_features that are positive
                    debug.best_cond = Inf;
                    for len=round(params.min_features/q):round(params.max_features/q)
                        G_rows = unique(sorted_rows(1:len, :));
                        G_sub = G(G_rows, :);
                        curr_cond = cond(G_sub);
                        Kon(len - round(params.min_features/q) + 1) = curr_cond;
                        if(curr_cond < debug.best_cond)
                            debug.best_cond = curr_cond;
                            debug.G = len;
                            final_features = G_rows;
                        end
                    end                    
                    
                case 'neq'
                    [~, sorted_rows] = sort(max(Stats, [], 2), 'descend');
                    
                    debug.best_cond = Inf;
                    for len=params.min_features:q:params.max_features                                                
                        H_sub = H(sorted_rows(1:len), :);
                        curr_cond = cond(H_sub);
                        if(curr_cond < debug.best_cond)
                            debug.best_cond = curr_cond;
                            debug.G = len;
                            final_features = sorted_rows(1:len);
                        end
                    end  
                    
                otherwise
                    error('Unknown selection method %s\n_sub', params.selection_method);
            end

        otherwise
            error('Unknown ranking method %s\n_sub', params.ranking_method);
    end
    
    final_features = selected_features(final_features);
    
%     selected_features = [];
%     tic;
%     switch params.method
%         case 'ttest'
%             [Classes,~ , pure_col_classes] = unique(annotations.pure.Class);
% 
%             for class = 1:numel(Classes)
%                 class_cols = find(pure_col_classes == class);
%                 null_cols = setdiff(1:numel(annotations.pure.Class), class_cols)';
%                 [pvalues] = mattest(H(:, class_cols), H(:, null_cols), 'Vartype', 'unequal');
%                 [~, q] = mafdr(pvalues);
%                 selected_rows = find(q <= qval_threshold & (mean(H(:, class_cols), 2) > mean(H(:, null_cols), 2)));
% %                 selected_genes = setdiff(unique(annotations.features.Primary_ID(selected_rows)), {''});
% %                 system(sprintf('mkdir -p output/sigGenes/%s', datasetName));
% %                 dlmcell(sprintf('output/sigGenes/%s/%s_qval=%.2e.txt', datasetName, Classes{class}, qval_threshold), selected_genes);
%                 selected_features = union(selected_features, selected_rows);
%             end            
% 
%         case 'ttest_cond'
%             [Classes,~ , pure_col_classes] = unique(annotations.pure.Class);
%             Fold = zeros(size(H, 1), numel(Classes));
%             selected_features = zeros(size(H, 1), 1);
%             total_selected = 0;
%             Class_mean = zeros(size(H, 1), numel(Classes));
%             Class_var = zeros(size(H, 1), numel(Classes));
%             for class = 1:numel(Classes)
%                 class_cols = find(pure_col_classes == class);
%                 Class_mean(:, class) = mean(H(:, class_cols), 2);
%                 Class_var(:, class) = var(H(:, class_cols), 0, 2);
%                 
%                 null_cols = setdiff(1:numel(annotations.pure.Class), class_cols)';
%                 [pvalues] = mattest(H(:, class_cols), H(:, null_cols), 'Vartype', 'unequal');
%                 [~, q] = mafdr(pvalues);
%                 fold_change = ( (mean(H(:, class_cols), 2) ./ mean(H(:, null_cols), 2)));
%                 selected_rows = find(q <= params.qval_threshold & params.fold_threshold < fold_change);
%                 
%                 Fold(total_selected+1:total_selected+numel(selected_rows), class) = fold_change(selected_rows);
%                 selected_features(total_selected+1:total_selected+numel(selected_rows)) = selected_rows;
%                 total_selected = total_selected + numel(selected_rows);


%                 
% %                 selected_genes = setdiff(unique(annotations.features.Primary_ID(selected_rows)), {''});
% %                 system(sprintf('mkdir -p output/sigGenes/%s', datasetName));
% %                 dlmcell(sprintf('output/sigGenes/%s/%s_qval=%.2e.txt', datasetName, Classes{class}, qval_threshold), selected_genes);
%             end                
%             Fold(total_selected+1:end, :) = [];
%             selected_features(total_selected+1:end) = [];
%             Class_mean_sub = Class_mean(selected_features, :); Class_mean_sub(total_selected+1:end, :) = [];
%             Class_var_sub = Class_var(selected_features, :); Class_var_sub(total_selected+1:end, :) = [];
%             
%             [min_marker_no, min_idx] = min(sum(logical(Fold)));
%             params.max_G = min(params.max_G, min_marker_no);
%             if(params.max_G < params.min_G)
%                 error('At least for celltype %s there are only %d markers which is smaller than G_min = %d. Try adjusting the qval_threshold to a larger value\n_sub', Classes{min_idx}, min_marker_no, params.min_G);
%             end
%             
%             [~, Fold_perm] = sort(Fold, 'descend');
%             sorted_rows = selected_features(Fold_perm);
%             
%             debug.best_cond = Inf;
%             for len=params.min_G:params.max_G
%                 H_rows = unique(sorted_rows(1:len, :));
%                 H_sub = H(H_rows, :);
%                 curr_cond = cond(H_sub);
%                 if(curr_cond < debug.best_cond)
%                     debug.best_cond = curr_cond;
%                     debug.G = len;
%                     selected_features = H_rows;
%                 end
%             end
%     end

end

