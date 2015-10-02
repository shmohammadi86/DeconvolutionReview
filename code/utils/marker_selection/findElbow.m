function [ maxDist, idxOfBestPoint ] = findElbow( v )
%% findElbow: Finds the index of element in input vector v that has the maximum distance to the line connecting the first <-> last element of the vector. 
% Warning: The input vector has to be pre-sorted in either increasing or decreasing order

    if(~isvector(v))
        maxDist =Inf;
        idxOfBestPoint = -1;
        return;
    end
    if(size(v, 2) == 1) 
        v = v';
    end
    
    % smooth the vector first
    v = fastsmooth(v, 10, 2, 1);
    
    
    nPoints = length(v);
    allCoord = [1:nPoints;v]';
    firstPoint = allCoord(1,:);

    %# get vector between first and last point - this is the line
    lineVec = allCoord(end,:) - firstPoint;

    %# normalize the line vector
    lineVecN = lineVec / sqrt(sum(lineVec.^2));

    %# find the distance from each point to the line:
    %# vector between all points and first point
    vecFromFirst = bsxfun(@minus, allCoord, firstPoint);

    
    %# To calculate the distance to the line, we split vecFromFirst into two 
    %# components, one that is parallel to the line and one that is perpendicular 
    %# Then, we take the norm of the part that is perpendicular to the line and 
    %# get the distance.
    %# We find the vector parallel to the line by projecting vecFromFirst onto 
    %# the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
    %# We project vecFromFirst by taking the scalar product of the vector with 
    %# the unit vector that points in the direction of the line (this gives us 
    %# the length of the projection of vecFromFirst onto the line). If we 
    %# multiply the scalar product by the unit vector, we have vecFromFirstParallel
    scalarProduct = dot(vecFromFirst, repmat(lineVecN,nPoints,1), 2);
    vecFromFirstParallel = scalarProduct * lineVecN;
    vecToLine = vecFromFirst - vecFromFirstParallel;

    
    %# distance to line is the norm of vecToLine
    distToLine = sqrt(sum(vecToLine.^2,2));


    %# now all you need is to find the maximum
    [maxDist, idxOfBestPoint] = max(distToLine);

end

