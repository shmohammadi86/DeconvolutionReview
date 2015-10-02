function [table, colNames] = my_dlmread( path, dlm)   
% Reads a tab-separated table.
% Input) path: input path the file to be read 
% Output) colNames: name of each column, parsed from the first row of the file, table: Rest of the lines are parsed as a string table.   
    if (nargin == 1)
        dlm = '\t';
    end
    
    fd = fopen(path, 'r');
    if(fd == -1)
        error('my_tdfread:: Cannot open file %s', path);
    else        
        header = fscanf(fd, '%[^\n]\n', 1);
        cols = textscan(header, '%s', 'Whitespace', '', 'Delimiter', dlm);
        colNames = cols{1};
        
        %single_format = sprintf('%%[^%s]', dlm);
        lin_table = textscan(fd,'%s', 'Delimiter', dlm, 'Whitespace', '');
        table = reshape(lin_table{1}, numel(colNames), numel(lin_table{1})/numel(colNames))';
    end
    fclose(fd);
end

