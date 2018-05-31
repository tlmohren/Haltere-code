%-------------------------------
% TMohren
% load comsol data and plot resulting strains
% 2017-08-09
%------------------------------
clc; clear all; close all
% addpathFolderStructureHaltere()
% 
sim.name = 'Haltere_CraneFly_Sphere';
% load(['data' filesep sim.name filesep sim.name '_allData'])
% 
% % Simulation parameters
% circleDistance = 300;               % distance from base to haltere 
% circleRadius = 150;                 % radius of haltere
% strainPoints = [circleDistance, 0,   circleRadius; ...
%                 circleDistance, 0,   - circleRadius; ...
%                circleDistance,  circleRadius, 0;...
%                circleDistance, -circleRadius, 0];
% 
% %% Analyze strain
% circleIndices = [];
% pointIndices = [];
% 
% % Find point and circle Indices
% pointIndices= findPointIndices( round(sim.strainXYZ,7) , strainPoints );
% circleIndices= findCircleIndices( round(sim.strainXYZ,7) , circleDistance,circleRadius);
% 
% %% sort circle Indices
% xyz = sim.strainXYZ;
% angle = atan2( xyz(3,circleIndices), xyz(2,circleIndices));
% angleDeg = rad2deg(angle)-180;
% angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
% [V,I_sort] = sort(angleDeg,'ascend');
% 
% Ind = circleIndices(I_sort);
% sidePoints = Ind( find( mod(V,90) == 0));
% 
% 
% %% 
% dotStyle = {'Marker','.'};
% % markerStyle = {'Marker','o','MarkerFaceColor','red','CData',[1,0,0]};
% % 
% figure()
%     plotMesh( sim.strainXYZ ,dotStyle) 
% %     plotMesh( sim.strainXYZ(:,circleIndices) ,markerStyle) 

model = createpde;
importGeometry(model,['data' filesep sim.name filesep 'meshText.stl'] );
pdegplot(model,'FaceLabels','on')

