function [informative_rows markers marker_scores] = getMostInformativeFeatures(A, varargin)
    function partitionMarkers(informative_rows)
            data = A(informative_rows, :);    
            [~, gene_row_idx] = sort(data, 2, 'descend');
            if(par.allowDuplicates)
                counts = round(sum(data, 2).^2 ./ sum(data.^2, 2));
                counts(isnan(counts)) = 0;
            else
                counts = ones(numel(data), 1);
            end
            
            marker_counts = zeros(celltype_no, 1);
            row_count = numel(informative_rows);
            for l = 1:celltype_no
                markers{l} = zeros(row_count, 1);
                marker_scores{l} = zeros(row_count, 1);
            end 
            
            for i = 1:numel(informative_rows)
                for j = 1:counts(i)             
                    marker_counts(gene_row_idx(i, j)) = marker_counts(gene_row_idx(i, j)) + 1;
                    markers{gene_row_idx(i, j)}(marker_counts(gene_row_idx(i, j))) = informative_rows(i);
                    marker_scores{gene_row_idx(i, j)}(marker_counts(gene_row_idx(i, j))) = overall_selectivity(informative_rows(i));
                end
            end   
            for l = 1:celltype_no
                markers{l}(marker_counts(l)+1:end) = [];
                marker_scores{l}(marker_counts(l)+1:end) = [];
            end             
    end

    params = inputParser;    
    
    params.addParamValue('sampling', 'uniform'      ,@(x) ischar(x)); % uniform/pooled
    params.addParamValue('selection_method', 'participation_ratio_stringent'      ,@(x) ischar(x));
    params.addParamValue('cut_threshold'           , 10      ,@(x) isscalar(x) & x <= 100 & x > 0); % Percentage/threshold to select top-ranked scores
    params.addParamValue('allowDuplicates', true     ,@(x) islogical(x)); 
    params.addParamValue('visualize', true     ,@(x) islogical(x)); 
    
    params.parse(varargin{:});
    
    par = params.Results;
    celltype_no = size(A, 2);
    gene_no = size(A, 1);

    % Compute tissue-selectivity of probes
    P = normalize(A, 'dim', 2, 'pnorm', 1); % Relative expression of each gene in different tissues
    I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene   
    H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
    overall_selectivity  = H ./ log2(celltype_no);    
    [~, overall_selectivity_idx] = sort(overall_selectivity);
    
    
    if(strcmp(par.selection_method, 'cond')) % Needs fixing
        min_idx = min(size(A));
        max_idx = min(10000, max(size(A))); % Limit to a max of 10K markers
        cuts = round( min_idx:10:max_idx );
        conds = smooth(arrayfun(@(x) cond(A(overall_selectivity_idx(1:x), :)), cuts), 25);        
        [~, min_idx] = min(conds);
        sigSelective_no = cuts(min_idx);
        par.sampling = 'uniform'; % Cond # method is only applicable to overall selectivity, not categorical.
    else        
        sigSelective_no = cut(overall_selectivity, 'selection_method', par.selection_method, 'selection_half', 'bottom', 'cut_threshold', par.cut_threshold);
    end

    Q = bsxfun(@plus, I,  H); % conditional specificity of the relative expression -- categorical  tissue-specificity of genes. . It has unit of "bits" and is between [0-2log2(k)]. Most genes with 0 ≤ Q ≤ 7 are tissue-specific
    categorical_selectivity = Q ./ (2*log2(celltype_no));

    markers = cell(celltype_no, 1);
    marker_scores = cell(celltype_no, 1);
    
    switch(par.sampling)
        case 'uniform'
            num_selected_rows = sigSelective_no;            
            informative_rows = overall_selectivity_idx(1:num_selected_rows);
            partitionMarkers(informative_rows);

        case 'pooled'
            num_selected_rows = min(gene_no, 2*sigSelective_no);                        
            informative_rows = overall_selectivity_idx(1:num_selected_rows);
            partitionMarkers(informative_rows);
            good_score = median(overall_selectivity(overall_selectivity_idx(1:sigSelective_no)));
            for k = 1:celltype_no
                celltype_selectivity_score = marker_scores{k};
                num_selected_celltype_markers = cut([good_score; celltype_selectivity_score], 'selection_method', par.selection_method, 'cut_threshold', par.cut_threshold)-1;            
                markers{k}(num_selected_celltype_markers+1:end) = [];
                marker_scores{k}(num_selected_celltype_markers+1:end) = [];
            end                       

        case 'categorical'
            for k = 1:celltype_no
                [scores, score_idx] = sort(categorical_selectivity(1:round(gene_no/2), k));
                num_selected_celltype_markers = cut(scores, 'selection_method', par.selection_method, 'selection_half', 'bottom');                    
                markers{k} = score_idx(1:num_selected_celltype_markers);
                marker_scores{k} = scores(1:num_selected_celltype_markers);      
            end
        otherwise
            error('getMostInformativeFeatures: Unknown sampling metod: %s', par.sampling);
    end
    
    informative_rows = [];
    for k = 1:celltype_no
        informative_rows = union(informative_rows, markers{k});
    end                           
end
