function [informative_rows, classes, class_IR, class_markers, class_markerNames, subclasses, subclass_IR, subclass_markers, subclass_markerNames] = getMarkers_multiResolution( A,  row_labels, column_labels, column_classes, varargin)
% TODO: 
% 1) Should we also use absolute value of expression in identifying marer genes?
% 2) Should we use quantile normlization within each class first? 
% 3) Should we use somthing rather than average?
% Input) 
% A: Expression matrix of reference cell types
% row_labels: Vector of unique id(s) corresponding to each row in A (probeset ID, gene name, ...)
% column_labels: Vector of unique id(s) corresponding to each column in A (Sample ID, ...)
% column_classes: A vector/cell array of class labels for each column in A (class, dup, ...)
% Output)
% informative_rows: Union of all within/between class informative rows
% Between class:
%        classes: Unique set of column classes in A
%        class_IR: Informative rows for betweeen class separation
%        class_markers: Classification of informative rows into
%        groups that are informative in different unique classes
%        class_markerNames: Row names (gene name, probeset Id, ...)
%        of informative markers in each unique class.
% Within class:
%        subclasses: Unique set of subclasses in each column class
%        subclass_IR: Informative rows for within subclass separation
%        subclass_markers: Classification of informative rows for the subclass into
%        groups that are informative in different subclasses
%        subclass_markerNames: Row names (gene name, probeset Id, ...)
%        of informative markers in each unique subclass.  
    params = inputParser;    
    
    params.addParamValue('sampling', 'uniform'      ,@(x) ischar(x)); % uniform/pooled
    params.addParamValue('selection_method', 'participation_ratio_stringent'      ,@(x) ischar(x));
    params.addParamValue('cut_threshold'           , 10      ,@(x) isscalar(x) & x <= 100 & x > 0); % Percentage/threshold to select top-ranked scores
    params.addParamValue('allowDuplicates', true     ,@(x) islogical(x)); 
    params.addParamValue('visualize', false     ,@(x) islogical(x)); 
    params.addParamValue('pnorm'         , 1     ,@(x) isscalar(x) & x > 0); 
    
    params.parse(varargin{:});
    par = params.Results;

    % Make sure all annotations are entered in vector form
    if(~isvector(column_labels) || ~isvector(row_labels))
        error('both column_labels and row_labels have to be a vector (either string cell or vector of integers)');           
    elseif(~isvector(column_classes))
        error('column_classes has to be a vector (either string cell or vector of integers)');
    end
    
    % Ensure both column vectors are oriented similarly
    if(size(column_classes, 2) == 1)
        column_classes =column_classes';
    end
    if(size(column_labels, 2) == 1)
        column_labels =column_labels';
    end
    
    % Ensure column labels and column groupings both math to ALL columns
    if(size(column_classes, 2) ~= size(A, 2) || size(column_labels, 2) ~= size(A, 2))
        error('column_classes and column_labels should have similar number of entities as the columns in the input matrix (A)');        
    end
    
    % If input classes are in form of an integer vector, convert it to cell
    % array of string
    if(~iscell(column_classes))
        column_classes = cellstr(num2str(column_classes'))';
    end
    
%     if(norm(sum(A)- ones(1, size(A, 2))) > 10^-6) % Make sure A is column normalized (if it already is, doesn't change anything)        
%         if(max(A(:)) <= 16) % If we normalized the log-transformed A we will NOT be able to use this condition !!
%             warning('Expression matrix A looks log-transformed. Please make sure to only use raw expression matrix in this function.');
%         end            
%         A = normalize(A, 'dim', 1, 'pnorm', par.pnorm); 
%     end
    
    [classes,~,ic] = unique(column_classes);
    class_no = numel(classes);
    
    subclasses=cell(class_no, 1); 
    subclass_IR=cell(class_no, 1); 
    subclass_markers=cell(class_no, 1); 
    subclass_markerNames=cell(class_no, 1); 
    A_averaged = zeros(size(A, 1), class_no);
    informative_rows = [];
    
    for class_id = 1:class_no
        col_idx = find(ic == class_id);
        A_averaged(:, class_id) = mean(A(:, col_idx), 2);
        subclasses{class_id} = column_labels(col_idx);
        if(1 < numel(col_idx))
            A_class = A(:, col_idx);            
            [subclass_IR{class_id}  subclass_markers{class_id}] = getMostInformativeFeatures(A_class, 'sampling', par.sampling, 'allowDuplicates', par.allowDuplicates, 'selection_method', par.selection_method, 'visualize', false);        
            subclass_markerNames{class_id} = cellfun(@(idx) row_labels(idx), subclass_markers{class_id}, 'UniformOutput', false);            
            informative_rows = union(informative_rows, subclass_IR{class_id});            
        end
    end

    [class_IR  class_markers] = getMostInformativeFeatures(A_averaged, 'sampling', par.sampling, 'allowDuplicates', par.allowDuplicates, 'selection_method', par.selection_method, 'visualize', par.visualize);        
    class_markerNames = cellfun(@(idx) row_labels(idx), class_markers, 'UniformOutput', false);                
    informative_rows = union(informative_rows, class_IR);                
end

