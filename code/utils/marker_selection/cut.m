function [ num_selected_entries ] = cut( v, varargin )
    params = inputParser;
    params.addParamValue('selection_method'           , 'participation_ratio',@(x) ischar(x) ); % selection method. Options:
    params.addParamValue('selection_half'           , 'top',@(x) ischar(x) ); % Should we select the top half or the bottom half of scores?
    params.addParamValue('cut_threshold'           , 100      ,@(x) isscalar(x) & x <= 100 & x > 0); % Percentage of top-ranked scores to choose    
        
    params.parse(varargin{:});
    par = params.Results;

    badV = isnan(v) | isinf(v);
    v(badV) = [];
    v = sort(v);    
    n = numel(v);
    
    switch(par.selection_method)
        case 'Elbow'
            [ ~, num_selected_entries ] = findElbow( v );            
            if(  (num_selected_entries < n/2 && strcmp(par.selection_half, 'top') ) || ... % Right-elbow 
                (num_selected_entries > n/2 && strcmp(par.selection_half, 'bottom') )  ) % Left-elbow
                    num_selected_entries = n - num_selected_entries;
            end
        case 'participation_ratio'
            num_selected_entries = round(sum(v)^2/sum(v.^2)); % participation ratio relaxed (allows more nonzeros)
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = n - num_selected_entries; % Participation ratio returns #nnz=#elements in Top half. If we want bottom half, we have to reverse it.
            end
        case 'participation_ratio_stringent'
            num_selected_entries = round(sum(v.^2)^2/sum(v.^4)); % participation ratio 
            if( strcmp(par.selection_half, 'bottom') )
                num_selected_entries = n - num_selected_entries; % Participation ratio returns #nnz=#elements in Top half. If we want bottom half, we have to reverse it.
            end
        case 'MeanDiff'
            q = nan(length(v)-1, 1);
            for i = 1:length(v)-1
                mu1 = mean(v(1:i)); mu2 = mean(v(i+1:end));
                sigma1 = std(v(1:i)); sigma2 = std(v(i+1:end));
                q(i) = (mu2 - mu1) / (sigma1 + sigma2);
            end
            [~, num_selected_entries] = min(q);       
        case 'AIC' % Unstable: needs work
            num_selected_entries = AIC_outlierDet(v, 0.9);
        case 'fix_threshold' 
            num_selected_entries = nnz(v >= par.cut_threshold);
        case 'fix_perc'
            num_selected_entries = round(par.cut_threshold*numel(v)/100);
        otherwise 
            error('Cell type selection method unknown: %s\n', par.selection_method);
    end  
end

