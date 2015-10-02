function [ M, G, annotations, known_proportions, H ] = loadDataset(paths)
%% [ M, G, annotations, known_proportions, H ] = loadDataset(paths); 
% Input) paths: an structure array constructed by setDatasetPaths()
% Output) M: mixture matrix, G: Reference profile matrix, annotation: structure
% array containing annotations for samples, references, features and pure
% samples, H: pure profile matrix
   %% Read M and G 
    fprintf('Loading mixture matrix (M) from %s ...\n', paths.expression.mixture); 
    [M, sample_id, feature_id_M] = my_tblread(paths.expression.mixture); 
    NaN_count = nnz(isnan(M(:)));
    if( NaN_count ~= 0) %tabbing error
        error('M has %d NaN values after importing. Are you sure it is properly formated (missing tab maybe?)', NaN_count);
    end
    % Check to make sure data is not log-transormed
    if(max(M(:)) <= 20) 
        warning('run:loadDataset', 'Expression matrix M looks log-transformed. Trying to reverse it here.');
        M = 2.^M;
    end    
    
    fprintf('Loading reference signature matrix (G) from %s ... \n', paths.expression.signature); 
    [G, reference_id, feature_id_G] = my_tblread(paths.expression.signature);     
    NaN_count = nnz(isnan(G(:)));
    if( NaN_count ~= 0) 
        error('G has %d NaN values after importing. Are you sure it is properly formated (missing tab maybe?)', NaN_count);
    end
    % Check to make sure data is not log-transormed
    if(max(G(:)) <= 20) % Affy has resolution of 2^16
        warning('run:loadDataset', 'Expression matrix G looks log-transformed. Trying to reverse it here.');
        G = 2.^G;        
    end
    
    common_features = intersect(feature_id_M, feature_id_G);
    [~, M_row_idx] = ismember(common_features, feature_id_M);
    M = M(M_row_idx, :);
    [~, G_row_idx] = ismember(common_features, feature_id_G);
    G = G(G_row_idx, :);
    
            
    %% Load annotations
        fprintf('Loading reference annotations from %s ...\n', paths.annotations.references);
        if(isempty(paths.annotations.references))
            fprintf('Missing reference annotations ... using first row of the signature expression file\n');
            annotations.references = table(reference_id', reference_id', reference_id', 'VariableNames', {'Reference_ID', 'Subclass', 'Class'}); % Assumes no grouping among references
        else
            annotations.references = readtable(paths.annotations.references, 'Delimiter', '\t');
        end

        fprintf('Loading sample annotations from %s ...\n', paths.annotations.samples);
        if(isempty(paths.annotations.samples))
            fprintf('Missing sample annotations ... using first row of the mixture expression file\n');
            annotations.samples= table(sample_id', 'VariableNames', {'Sample_ID'});
        else
            annotations.samples =  readtable(paths.annotations.samples, 'Delimiter', '\t');
        end

        fprintf('Loading feature annotations from %s ...\n', paths.annotations.features);
        if(isempty(paths.annotations.features))
            fprintf('Missing feature annotations ... using first column of the mixture expression file\n');
            annotations.features= table(common_features, common_features, 'VariableNames', {'Feature_ID', 'Gene_Symbol'});
        else            
            annotations.features = readtable(paths.annotations.features, 'Delimiter', '\t');
            Max_Match_count = 0; Max_Match_idx = 0;

            for i = 1:size(annotations.features, 2)
                match_count = nnz(ismember(common_features, table2cell(annotations.features(:, i))));
                if(Max_Match_count < match_count )
                    Max_Match_count = match_count;
                    Max_Match_idx = i;
                    break;
                end
            end
                
            if(0 < Max_Match_count)
                base_table = table(common_features, 'VariableNames', annotations.features.Properties.VariableNames(Max_Match_idx));
                annotations.features = outerjoin(base_table, annotations.features,'MergeKeys',true,'Type','left');
                annotations.features.Properties.VariableNames{1} = 'Feature_ID';
            else
                annotations.features= table(common_features, common_features, 'VariableNames', {'Feature_ID', 'Gene_Symbol'});
            end                    
        end
        primary_cols = find(ismember(upper(annotations.features.Properties.VariableNames), upper({'GeneName', 'Gene_Name', 'GeneSymbol', 'Gene_Symbol'})));        
        if(isempty(primary_cols))
            error('Cannot find gene symbol annotations to set as primary feature ID');
        end
        annotations.features.Properties.VariableNames{primary_cols(1)} = 'Primary_ID';

    %% Loading pure profiles, if exists
    if(~isempty(paths.expression.pure)) % Is not mandatory, but, if provided, make sure the number of features match features in G
        fprintf('Loading pure signature matrix (H) from %s ... \n', paths.expression.pure); 
        [H, ~, feature_id_Z] = my_tblread(paths.expression.pure);     
        if( ~isempty(setdiff(feature_id_Z, feature_id_G)) || ~isempty(setdiff(feature_id_G, feature_id_Z)) )
            error('Features in G (signature profile) and H (pure profile) should match.');
        end
        NaN_count = nnz(isnan(H(:)));
        if( NaN_count ~= 0) 
            error('H has %d NaN values after importing. Are you sure it is properly formated (missing tab maybe?)', NaN_count);
        end
        % Check to make sure data is not log-transormed
        if(max(H(:)) <= 20) % Affy has resolution of 2^16
            warning('run:loadDataset', 'Expression matrix G looks log-transformed. Trying to reverse it here.');
            H = 2.^H; 
        end
        H = H(G_row_idx, :);
        
        annotations.pure =  readtable(paths.annotations.pure, 'Delimiter', '\t');        
    else
        H = [];
        annotations.pure= table();
    end
%% Summarize G and sort alphabetically based on the reference name
    reference_types = setdiff(unique(annotations.references.Type), '-'); % if more one reference/ cell type is provided, compute their average over unique values of references.Type
    if(numel(reference_types) ~= size(G ,2)) % summarize G
        new_L = zeros(size(G, 1), numel(reference_types));
        for i = 1:numel(reference_types)
            selected_cols = find(strcmp(reference_types{i}, annotations.references.Type));
            if(numel(selected_cols) == 1)
                new_L(:, i) = G(:, selected_cols);
            else
                new_L(:, i) = mean(G(:, selected_cols), 2);
            end
        end
        G = new_L;        
    else
        [~, col_perm] = ismember(reference_types, annotations.references.Type);
        G = G(:, col_perm);
    end

    [~, row_idx] = ismember(reference_types, annotations.references.Type);
    annotations.references = annotations.references(row_idx, :);    

    %% Load known coefficients fo controlled experiments
    if( ~isempty(paths.coeff) )        
        [known_proportions.C, ~, known_proportions.types]=my_tblread(paths.coeff);         
        if(size(known_proportions.C, 2) ~= size(M, 2))
            error('Number of samples in M (%d) does not match the given known percentages (%d)\n', size(Pk, 2), size(M, 2));
        end
        
        % Make sure columns of G and rows of C are in the same order
        [~, row_idx] = ismember(reference_types, known_proportions.types);
        known_proportions.C = known_proportions.C(row_idx, :);
        known_proportions.types = known_proportions.types(row_idx, :);
        
        known_proportions.C = normalize(known_proportions.C, 'pnorm', 1, 'dim', 1);
    end      

end

