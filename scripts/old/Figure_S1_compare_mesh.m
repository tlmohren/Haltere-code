clc;clear all;close all
addpathFolderStructureHaltere()

load(['data' filesep 'Cranefly_SphereVsEllipsoid_FEAresults'],'FEA')

%% Determine Circle locations
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   

xMatch = find( round( FEA(1).XYZ(:,1), 7) == circleDistance );
yMatch = find( round( abs( FEA(1).XYZ(:,2) ), 7) == circleRadius );
zMatch = find( round( abs( FEA(1).XYZ(:,3) ), 7) == circleRadius );

sideInds = intersect(xMatch,yMatch);
topInds = intersect(xMatch,zMatch);

%% 


figure();
for j = 1:9
    subplot(3,3,j)
    scatter3( FEA(j).XYZ(:,1), FEA(j).XYZ(:,2), FEA(j).XYZ(:,3)  )
    axis equal
    
end
