% plots effect of rotation on strain at haltere sides 
clc;clear all;close all
addpathFolderStructureHaltere()

FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';

for j =  1:length(FEA)
    tic
    [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
    toc 
end

% Determine Circle locations
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   

xMatch = find( FEA(1).xyz(:,1) == circleDistance );
yMatch = find( abs( FEA(1).xyz(:,2) ) == circleRadius );
zMatch = find( abs( FEA(1).xyz(:,3) ) == circleRadius );

sideInds = intersect(xMatch,yMatch);
topInds = intersect(xMatch,zMatch);

% Plot figure
t = 0:0.001:0.35;
figure()
subplot(211); hold on
    plot(t,FEA(1).strain(topInds,:))
    plot(t,FEA(2).strain(topInds,:))
    xlabel('Time (s)');
    ylabel('$\epsilon$')
    legend( '$\Omega = 0$ Top', '$\Omega = 0$ Bottom','$\Omega = 10$ Top' ,'$\Omega = 10$ Bottom', 'Location', 'NorthEastOutside')
subplot(212); hold on
    plot(t, FEA(1).strain(sideInds,:))
    plot(t, FEA(2).strain(sideInds,:))
    xlabel('Time (s)');
    ylabel('$\epsilon$')
    legend( '$\Omega = 0$ Left', '$\Omega = 0$ Right','$\Omega = 10$ Left' ,'$\Omega = 10$ Right', 'Location', 'NorthEastOutside')