function [ pointIndices ] = findPointIndices( xyz, strainPoints )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    for j = 1:size(strainPoints,1)     % for the 4 points around the haltere (top, bottom, left,right) 
        [~,yBin] = find( xyz(2,:) == strainPoints(j,2));
        [~,zBin] = find( xyz(3,:) == strainPoints(j,3));
        yzBin = intersect(yBin,zBin);   
        [~,yzBinMin] = min(abs( strainPoints(j,1)  - xyz(1,yzBin) ) );
        min(abs( strainPoints(j,1)  - xyz(1,yzBin) ) );
        pointIndices(j) = yzBin(yzBinMin);                              % loc_I is used for base points 
    end

end

