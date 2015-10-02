function [table, colNames, rowNames] = my_tblread2( path, dlm)
    if (nargin == 1)
        dlm = '\t';
    end
    
    fd = fopen(path, 'r');
    if(fd == -1)
        error('my_tblread:: Cannot open file %s', path);
    else        
        %header = textscan(fd, '%[^\n]\n', 1);
        header = fgets(fd);
        cols = textscan(header, '%s', 'Whitespace', '', 'Delimiter', dlm);
        colNames = cols{1}(1:end)';

        FormatString = sprintf('%%s%s', repmat('%f',1, numel(colNames)));
        T = textscan(fd, FormatString, 'CollectOutput', 1, 'Whitespace', '', 'Delimiter', dlm);
        
        rowNames = T{1};
        table = T{2};
    end
    fclose(fd);
end

