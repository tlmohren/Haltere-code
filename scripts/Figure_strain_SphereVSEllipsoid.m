% plots effect of rotation on strain at haltere sides 
clc;clear all;close all
addpathFolderStructureHaltere()

FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
FEA(3).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
FEA(4).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
FEA(5).name = 'Haltere_CraneFly_ellipsoidVer_Om0';
FEA(6).name = 'Haltere_CraneFly_ellipsoidVer_Om10';

for j =  1:length(FEA)
    tic
    [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
    toc 
end
%% 
% Determine Circle locations
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   

xMatch = find( FEA(1).xyz(:,1) == circleDistance );
yMatch = find( abs( FEA(1).xyz(:,2) ) == circleRadius );
zMatch = find( abs( FEA(1).xyz(:,3) ) == circleRadius );

sideInds = intersect(xMatch,yMatch);
topInds = intersect(xMatch,zMatch);

%% Plot stuff 
t = 0:0.001:0.35;
figure()
for j = 1:3
subplot(3,2,j*2-1); hold on
    plot(t,FEA(j*2-1).strain(topInds,:))
    plot(t,FEA(j*2).strain(topInds,:))
    xlabel('Time (s)');
    ylabel('$\epsilon$')
    
subplot(3,2,j*2); hold on
    plot(t, FEA(j*2-1).strain(sideInds,:))
    plot(t, FEA(j*2).strain(sideInds,:))
    xlabel('Time (s)');
    ylabel('$\epsilon$')
end

%% 
figure()
for j = 1:3
subplot(1,2,1); hold on
    plot(t,FEA(j*2-1).strain(topInds,:))
    plot(t,FEA(j*2).strain(topInds,:))
    xlabel('Time (s)');
    ylabel('$\epsilon$')
    
subplot(1,2,2); hold on
    plot(t, FEA(j*2-1).strain(sideInds,:))
    plot(t, FEA(j*2).strain(sideInds,:))
    xlabel('Time (s)');
    ylabel('$\epsilon$')
end