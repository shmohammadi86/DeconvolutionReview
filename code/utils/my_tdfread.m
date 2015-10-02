function [table colNames] = my_tdfread( path)   
% Reads a tab-separated table.
% Input) path: input path the file to be read 
% Output) colNames: name of each column, parsed from the first row of the file, table: Rest of the lines are parsed as a string table.   
    fd = fopen(path, 'r');
    if(fd == -1)
        error('my_tdfread:: Cannot open file %s', path);
    else        
        header = fscanf(fd, '%[^\r\n]\r\n', 1);
        cols = textscan(header, '%s', 'Whitespace', '', 'Delimiter', '\t');
        colNames = cols{1};

        FormatString = sprintf('%s%%[^\\r\\n]', repmat('%[^\t]\t',1, numel(colNames)-1));
        table = textscan(fd, FormatString, 'Whitespace', '', 'Delimiter', '\t');
    end
    fclose(fd);
end

