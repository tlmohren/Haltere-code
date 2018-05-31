clc;clear all;close all

addpathFolderStructureHaltere()
tic 
% load(['data' filesep 'collected_strainYYData'],'sim')
% load(['data' filesep 'collected_longRun'],'sim')
load(['data' filesep 'Cranefly_Sphere_fling_strainDeform'],'sim')

toc
%% plot mesh

sim(1).deform = sim(1).deform/max(sim(1).deform(:))*10;

mode = 2;
figure()
    plotStyle = {'Marker','.'};
    cData = sim(1).strain(:,mode);
    plotMeshColor( sim(1).XYZ, cData, plotStyle );
    colorbar
%     caxis([-1,1])
    
    
    
figure()
    plotStyle = {'Marker','.'};
    cData = sim(1).deform(:,mode,1);
    plotMeshColor( sim(1).XYZ, cData, plotStyle );
    colorbar
    
%     
% figure()
%     plotStyle = {'Marker','.'};
%     cData = sim(1).deform(:,mode,2);
%     plotMeshColor( sim(1).XYZ, cData, plotStyle );
%     colorbar
% figure()
%     plotStyle = {'Marker','.'};
%     cData = sim(1).deform(:,mode,3);
%     plotMeshColor( sim(1).XYZ, cData, plotStyle );
%     colorbar
%     caxis([-1,1])
%% Find points of interest 
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   

circleIndices= findCircleIndices( round(sim(1).XYZ,7) , circleDistance,circleRadius);
% xyz = sim(1).XYZ;
angle = atan2( sim(1).XYZ(3,circleIndices), sim(1).XYZ(2,circleIndices));
angleDeg = rad2deg(angle)-90;
angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
[V,I_sort] = sort(angleDeg,'ascend');
Ind = circleIndices(I_sort);
sidePoints = Ind( find( mod( round(V,3),90) == 0));




edgeIndices(:,1)= find( round(sim(1).XYZ(2,:),7) == circleRadius);
% edgeIndices(:,2)= find( round(sim(1).XYZ(2,:),7) == -circleRadius);
edgeIndices(:,2)= find( round(sim(1).XYZ(2,:),7) == -circleRadius);
edgeIndices(:,3)= find( round(sim(1).XYZ(3,:),7) == circleRadius);
edgeIndices(:,4)= find( round(sim(1).XYZ(3,:),7) == -circleRadius);
% xyz = sim(1).XYZ;
angle = atan2( sim(1).XYZ(3,circleIndices), sim(1).XYZ(2,circleIndices));
%% 
mode = 5


x =  sim(1).XYZ(1,edgeIndices(:,1) );
y = sim(1).deform( edgeIndices(:,1),mode,2);
z = sim(1).deform( edgeIndices(:,2),mode,3);



figure();
subplot(211);hold on

plot( sim(1).XYZ(1,edgeIndices(:,1)  ), sim(1).deform( edgeIndices(:,1),mode,3) )
plot( sim(1).XYZ(1,edgeIndices(:,1)  ), sim(1).deform( edgeIndices(:,2),mode,3) ,'o')
subplot(212); hold on
plot( sim(1).XYZ(1,edgeIndices(:,2)  ), sim(1).deform( edgeIndices(:,3),mode,3) )
plot( sim(1).XYZ(1,edgeIndices(:,2)  ), sim(1).deform( edgeIndices(:,4),mode,3) ,'o' )

% dotStyle = {'Marker','.'};
% markerStyle = {'Marker','o','MarkerFaceColor','red','CData',[1,0,0]};
% % figure()
%     plotMesh( sim(1).XYZ ,dotStyle) 
%     plotMesh( sim(1).XYZ (:,sidePoints) ,markerStyle)
%     
%     xlabel('X'); ylabel('Y'); zlabel('Z')
%     
    %% plot strains
% dt = 0.001;
% 
% lineStyle  = {'-','-'};
% figure()
% for j = 1:length(sim)
%     t = 0: dt: (length(sim(j).phi)-1)*dt;
%     subplot(211); hold on
%     plot( t, sim(j).strain(:,sidePoints([1,3]))  ,lineStyle{j})
%     subplot(212); hold on 
%     plot( t, sim(j).strain(:,sidePoints([2,4]))  ,lineStyle{j})
% end
% % plot( t, sim(2).strain(:,sidePoints,4)  ,'o-')
% subplot(211); xlabel('Time (s)'); ylabel('strain in x')
% subplot(212); xlabel('Time (s)'); ylabel('strain in x')



%% deformation analysis 

% find points of interest 
circleDistance = 4580;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   

j = 1;
circleIndices= findCircleIndices( round(sim(j).XYZ,7) , circleDistance,circleRadius);
% xyz = sim(1).XYZ;
angle = atan2( sim(j).XYZ(3,circleIndices), sim(j).XYZ(2,circleIndices));
angleDeg = rad2deg(angle)-90;
angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
[V,I_sort] = sort(angleDeg,'ascend');
Ind = circleIndices(I_sort);
% sidePoints = Ind( find( mod(V,90) == 0));
sidePoints = Ind( find( mod( round(V,3),90) == 0));

dotStyle = {'Marker','.','MarkerFaceColor','k'};
markerStyle = {'Marker','o','MarkerFaceColor','red','CData',[1,0,0]};
% figure()
%     plotMesh( sim(j).XYZ ,dotStyle) 
%     plotMesh( sim(j).XYZ (:,sidePoints) ,markerStyle)
    
%     sidePoints
%% 
% figure(); hold on
% for j = 1:3
%     subplot(3,1,j); hold on
%     plot( sim(1).deform(:, sidePoints,j) )
% %     plot( sim(2).deform(:, sidePoints,j) )
% end

%%

[n_points, n_eigs, n_disps ] = size(sim(1).deform);
for j = 1:length(sim)
    for k = 1:n_eigs
    tic
    % get angles / size
%     phi = sim(j).phi;
%     theta = -sim(j).theta; 
    [ n_points, n_eigs, n_deform] = size(sim(j).deform);
    
%     % compute actual point locations 
    sim(j).xyzPoints = sim(j).deform(:,:,2:end) + ...       
        permute( repmat(sim(j).XYZ(:,:)',1,1, n_eigs), [1,3,2] ) ;

    % transpose points into rigid haltere frame 

%     sim(j).xyzHaltereFrame = sim(j).deform
    % compute middle point 
    sim(j).xyzHaltereFrameMiddle = squeeze(  mean( sim(j).xyzPoints ( sidePoints([1,3]),:), 2)  );

    % dz = ( xyzHaltereFrame(:,sidePoints([2]),3) - xyzHaltereFrame(:,sidePoints([4]),3) )'
    sim(j).dz(k) = diff( sim(j).xyzPoints ( sidePoints([2,4]),k,3)');
    sim(j).dy(k) = diff( sim(j).xyzPoints ( sidePoints([2,4]),k,2)');
    sim(j).dx(k) = diff( sim(j).xyzPoints ( sidePoints([2,4]),k,1)');
%     sim(j).dxMiddle = diff( sim(j).xyzHaltereFrameMiddle(:,1)');
%     sim(j).dzTop = diff( sim(j).xyzHaltereFrame(:,sidePoints([1,3]),3)');
%     sim(j).dyTop = diff( sim(j).xyzHaltereFrameMiddle(:,sidePoints([1,3]),2)');
%     sim(j).dyTop = diff( sim(j).xyzHaltereFrame(:,sidePoints([1,3]),1)');

%     sim(j).twistAngle(k) = atan2(sim(j).dz,sim(j).dy);
    sim(j).twistAngle(k) = atan2(sim(j).dz(k),sim(j).dy(k));
%     sim(j).yAngle = atan( sim(j).xyzHaltereFrameMiddle(:,3) ./ sim(j).xyzHaltereFrameMiddle(:,1));
%     sim(j).zAngle = atan( sim(j).xyzHaltereFrameMiddle(:,2) ./ sim(j).xyzHaltereFrameMiddle(:,1));
%     sim(j).twistAngleTop = atan( sim(j).dyTop ./ sim(j).dzTop );
    toc 
    end
end


%% 