clc;clear all;close all
addpathFolderStructureHaltere()

FEA(1).name = 'Haltere_CraneFly_Sphere_eigenfrequency';
FEA(2).name = 'Haltere_CraneFly_ellipsoidHor_eigenfrequency';
FEA(3).name = 'Haltere_CraneFly_ellipsoidVer_eigenfrequency';

parameters = { 'disp (µm)','u (µm)','v (µm)','w (µm)','freq (1/s)'}; % select which parameters to load 

for j =  1:length(FEA)
    tic
    [FEA(j).xyz, FEA(j).data, FEA(j).header] = loadCSV( ['data' filesep  FEA(j).name], parameters);
    toc 
end

%% Determine Circle locations

for j = 1:length(FEA)
    circleDistance = 4500;               % distance from base to haltere 
    circleRadius = 150;                 % radius of haltere   
    
    xMatch = find( round( abs( FEA(j).xyz(:,1) ), 3) == circleDistance );
    yMatch = find( round( abs( FEA(j).xyz(:,2) ), 3) == circleRadius );
    zMatch = find( round( abs( FEA(j).xyz(:,3) ), 3) == circleRadius );

    FEA(j).sideInds = intersect(xMatch,yMatch);
    FEA(j).topInds = intersect(xMatch,zMatch);
    FEA(j).eigenfrequencies = FEA(j).data(1,:,5);
end
%
j = 1
FEA(j)
figure()
        scatter3( FEA(j).xyz(:,1) , FEA(j).xyz(:,2), FEA(j).xyz(:,3), '.k')
        hold on
        scatter3( FEA(j).xyz(FEA(j).sideInds ,1) , FEA(j).xyz(FEA(j).sideInds ,2), FEA(j).xyz(FEA(j).sideInds ,3),100,'or','filled')
        scatter3( FEA(j).xyz(FEA(j).topInds ,1) , FEA(j).xyz(FEA(j).topInds ,2), FEA(j).xyz(FEA(j).topInds ,3),100,  'ob','filled')
    axis equal
    xlabel('x');ylabel('y');zlabel('z')
    
%% find eigenfrequencies

enlarged = 2e3;
for j = 1:length(FEA)
    deforms = FEA(j).data(:,:,[2,3,4]) / max(max(max(FEA(j).data(:,:,[2,3,4])))) ; 
    [ ~ ,n_times, ~] = size( FEA(j).data(:,:,[2,3,4]) );
    
    FEA(j).xyzPoints = deforms*enlarged  + ...       
            permute( repmat( FEA(j).xyz,1,1, n_times), [1,3,2] ) ;
        
    % compute deformations
    FEA(j).xyzHaltereFrameMiddle = squeeze(  mean( FEA(j).xyzPoints( FEA(j).sideInds,:,:), 1)  );
    FEA(j).dz = diff( FEA(j).xyzPoints( FEA(j).sideInds,:,3));
    FEA(j).dy = diff( FEA(j).xyzPoints( FEA(j).sideInds,:,2));
%     
    FEA(j).twistAngle = atan2(FEA(j).dz, abs( FEA(j).dy) );
    FEA(j).yAngle = atan( FEA(j).xyzHaltereFrameMiddle(:,3) ./ FEA(j).xyzHaltereFrameMiddle(:,1));
    FEA(j).zAngle = atan( FEA(j).xyzHaltereFrameMiddle(:,2) ./ FEA(j).xyzHaltereFrameMiddle(:,1));
end

%% 
for j = 1:length(FEA)
    figure()
    for k = 1:6
        subplot(3,2,k); hold on 
        scatter3( FEA(j).xyz(:,1) , FEA(j).xyz(:,2), FEA(j).xyz(:,3), '.k' ,'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
        scatter3( FEA(j).xyzPoints(:,k,1) , FEA(j).xyzPoints(:,k,2), FEA(j).xyzPoints(:,k,3),'.r')
        title(['Frequency = ',num2str( round( FEA(j).eigenfrequencies(k),0) ), newline ,...
            '$\delta \phi$ = ', num2str( round( FEA(j).yAngle(k),3) ),...
            ', $\delta \theta$ = ', num2str(round( FEA(j).zAngle(k) ,3) ),...
            ', $\delta \gamma$ = ', num2str(round( FEA(j).twistAngle(k) ,3) ),...
            ])
        axis equal
        view(46,21)
        axis([0,7000,-2e3,2e3,-2e3,2e3])
        axis off
    end
end
