function markers = getMarkers(A, informative_rows)
    data = A(informative_rows, :);    
    n = size(data, 1); % number of rows/probesets/features
    m = size(data, 2); % number of cols/references
    
    [~, idx] = sort(data, 2, 'descend');
    counts = round(sum(data, 2).^2 ./ sum(data.^2, 2));
    
    markers = cell(m, 1);
    for i = 1:numel(informative_rows)
        for j = 1:counts(i)             
            markers{idx(i, j)} = [markers{idx(i, j)}, informative_rows(i)];
        end
    end
    
end

