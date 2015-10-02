function [ dataset_looking_output ] = dummy_dataset( varargin)
    dataset_looking_output = struct();
    if(ischar(varargin{1}) && strcmp(varargin{1}, 'File')) 
        path = varargin{2};        
        fprintf('Dummy DS: Reading from file %s ... \n', path);
        [vars var_names] = my_tdfread(path); 
        for i = 1:numel(var_names)
            dataset_looking_output.(var_names{i}) = vars{i};
        end                
    else
        fprintf('Use given variable ... \n');
        var_names = varargin{end};
        for i = 1:numel(var_names)
            dataset_looking_output.(var_names{i}) = varargin{i};
        end
    end
end

