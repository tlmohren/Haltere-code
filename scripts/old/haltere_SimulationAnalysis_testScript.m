clc;clear all;close all

addpathFolderStructureHaltere()
tic 
% load(['data' filesep 'collected_strainYYData'],'sim')
% load(['data' filesep 'collected_longRun'],'sim')
load(['data' filesep 'Cranefly_SphereK0_00001compare'],'sim')

toc
%% plot mesh
figure()
    plotStyle = {'Marker','.'};
    cData = sim(1).strain(50,:);
    plotMeshColor( sim(1).XYZ, cData, plotStyle );
    colorbar

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



% dotStyle = {'Marker','.'};
% markerStyle = {'Marker','o','MarkerFaceColor','red','CData',[1,0,0]};
% % figure()
%     plotMesh( sim(1).XYZ ,dotStyle) 
%     plotMesh( sim(1).XYZ (:,sidePoints) ,markerStyle)
%     
%     xlabel('X'); ylabel('Y'); zlabel('Z')
%     
    %% plot strains
dt = 0.001;

lineStyle  = {'-','-'};
figure()
for j = 1:length(sim)
    t = 0: dt: (length(sim(j).phi)-1)*dt;
    subplot(211); hold on
    plot( t, sim(j).strain(:,sidePoints([1,3]))  ,lineStyle{j})
    subplot(212); hold on 
    plot( t, sim(j).strain(:,sidePoints([2,4]))  ,lineStyle{j})
end
% plot( t, sim(2).strain(:,sidePoints,4)  ,'o-')
subplot(211); xlabel('Time (s)'); ylabel('strain in x')
subplot(212); xlabel('Time (s)'); ylabel('strain in x')

%% freq analysis 
mult = 1;
len = 250;

lineStyle  = {'-','--'};
figure()
for j = 1:length(sim)
    subplot(211); hold on
    strainPat = sim(j).strain(end-len+1:end,sidePoints([1,3])); 
    sig = repmat(strainPat,10,1);
    strainPat = interp1( 1:len,strainPat,1:1/mult:len);
    sig = repmat(strainPat,mult,1);
%     sig = sim(j).strain(end-249:end,sidePoints([1,3])); 
    
%     [k,y] = quick_fft(sig,1/dt);
    [k,y] = quick_fft(sig,1/dt*mult);
    plot(k,y,lineStyle{j})
%     aa = find(y>max(y)*1e-9);
%     plot(k(aa),y(aa),'*-')
    
    subplot(212); hold on 
    strainPat = sim(j).strain(end-len+1:end,sidePoints([2,4])); 
    strainPat = interp1( 1:len,strainPat,1:1/mult:len);
    sig = repmat(strainPat,mult,1);
%     sig = sim(j).strain(end-249:end,sidePoints([2,4])); 
    [k,y] = quick_fft(sig,1/dt*mult);
    
    plot(k,y,lineStyle{j})
%     aa = find(y>max(y)*1e-9);
%     plot(k(aa),y(aa),'*-')
    
end
%     axOpts = {'ALim',{}};
% %     axOpts = {'XLIM',[0,250]}; 
%     subplot(211); ax = gca();set(ax,axOpts{:});
%     subplot(212);  ax = gca();set(ax,axOpts{:});

%% deformation analysis 

% find points of interest 
circleDistance = 4580;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   

j = 2;
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
for j = 1:length(sim)
    tic
    % get angles / size
    phi = sim(j).phi;
    theta = -sim(j).theta; 
    [n_times, n_points, n_deform] = size(sim(j).deform);
    
    % compute actual point locations 
    sim(j).xyzPoints = sim(j).deform(:,:,:) + ...       
        permute( repmat(sim(j).XYZ(:,:)',1,1, n_times), [3,1,2] ) ;

    % transpose points into rigid haltere frame 
    for k = 1: n_times
        for l = 1:n_points
            eul_theta = [cos(theta(k))   sin(theta(k))    0 ; ...
                         -sin(theta(k)) cos(theta(k))     0 ;...
                         0          0               1]^-1;
            eul_phi = [cos(phi(k))      0       sin(phi(k)) ; ...
                        -0              1       0 ;...
                        -sin(phi(k))    0        cos(phi(k))]^-1;
%             sim(j).xyzHaltereFrame(k,l,:) = eul_theta*eul_phi* squeeze( sim(j).xyzPoints(k,l,:) ) ;
            sim(j).xyzHaltereFrame(k,l,:) = eul_phi*eul_theta* squeeze( sim(j).xyzPoints(k,l,:) ) ;
        end
    end
    % compute middle point 
    sim(j).xyzHaltereFrameMiddle = squeeze(  mean( sim(j).xyzHaltereFrame(:, sidePoints([1,3]),:), 2)  );

    % dz = ( xyzHaltereFrame(:,sidePoints([2]),3) - xyzHaltereFrame(:,sidePoints([4]),3) )'
    sim(j).dz = diff( sim(j).xyzHaltereFrame(:,sidePoints([2,4]),3)');
    sim(j).dy = diff( sim(j).xyzHaltereFrame(:,sidePoints([2,4]),2)');
%     sim(j).dxMiddle = diff( sim(j).xyzHaltereFrameMiddle(:,1)');
%     sim(j).dzTop = diff( sim(j).xyzHaltereFrame(:,sidePoints([1,3]),3)');
%     sim(j).dyTop = diff( sim(j).xyzHaltereFrameMiddle(:,sidePoints([1,3]),2)');
%     sim(j).dyTop = diff( sim(j).xyzHaltereFrame(:,sidePoints([1,3]),1)');
    sim(j).twistAngle = atan2(sim(j).dz,sim(j).dy);
    sim(j).yAngle = atan( sim(j).xyzHaltereFrameMiddle(:,3) ./ sim(j).xyzHaltereFrameMiddle(:,1));
    sim(j).zAngle = atan( sim(j).xyzHaltereFrameMiddle(:,2) ./ sim(j).xyzHaltereFrameMiddle(:,1));
%     sim(j).twistAngleTop = atan( sim(j).dyTop ./ sim(j).dzTop );
    toc 
end
%% 

% thetaTest  = atan2( sim(2).xyzHaltereFrameMiddle(:,2) , sim(2).xyzHaltereFrameMiddle(:,1) );  
% 
%  figure();plot(-theta); hold on; plot(thetaTest)
%  legend('$\theta$','$\theta$ corrected')
%  

%% Look at deformations 
figure()
for j = 1:length(sim)
    subplot(421); hold on 
    plot(sim(j).xyzHaltereFrameMiddle(:,1))
    subplot(423);hold on 
    plot(sim(2).xyzHaltereFrameMiddle(:,2))
    subplot(425);hold on 
    plot(sim(j).xyzHaltereFrameMiddle(:,3))
    subplot(427);hold on 
    plot( sim(j).twistAngle)
%     plot( sim(j).twistAngleTop,'-')
end

dt = 0.001;
mult = 1;
fftLen = 200;
% rangeWords = 'end-fftLen+1:end'
% figure()
for j = 1:length(sim)
    for k = 1:3
        subplot(4,2,2*k); hold on
        sig = sim(j).xyzHaltereFrameMiddle(end-fftLen+1:end,k);
        sig = sim(j).xyzHaltereFrameMiddle(end-fftLen+1:end,k);
        [k,y] = quick_fft(sig,1/dt*mult);
        plot(k,y,'*-')
    end
        subplot(428); hold on
        sig =sim(j).twistAngle(end-fftLen+1:end);
        [k,y] = quick_fft(sig,1/dt*mult);
        plot(k,y,'*-')
        size(sig)
end















%% check deformations 

% figure()
% 
% for j = 1:length(sim)
%     subplot(211);hold on 
%     plot(  sim(j).xyzHaltereFrame(:,sidePoints([2,4]),3))
%     subplot(212);hold on 
% %     plot(   sim(1).xyzHaltereFrame(:,sidePoints([4]),3)  - sim(1).xyzHaltereFrame(:,sidePoints([2]),3)  )
%     plot(   diff( sim(j).xyzHaltereFrame(:,sidePoints([2,4]),3)' )  )
% 
% end
 
%% 
% 
% 
% figure()
% 
% for j = 1:length(sim)
%     subplot(211);hold on 
%     plot( sim(j).dy)
%     subplot(212);hold on 
%     plot( sim(j).dz)
% end
 

%% 
figure()
for j = 1:2
    subplot(3,2,1) ; hold on
    plot( sim(j).yAngle  )
    subplot(3,2,3)  ; hold on
    plot( sim(j).zAngle  )
    subplot(3,2,5)  ; hold on
    plot( sim(j).twistAngle   )
end
% subplot(321); ylabel('In-plane angle','Rotation',0)
% subplot(323); ylabel('Out-plane angle','Rotation',0)
% subplot(325); ylabel('Twist angle','Rotation',0)%% 

% figure()
for j = 1:2
    subplot(3,2,2)  ; hold on
    plot( sim(j).phi, sim(j).yAngle  )
    subplot(3,2,4)  ; hold on
    plot( sim(j).phi, sim(j).zAngle  )
    subplot(3,2,6)  ; hold on
    plot( sim(j).phi,sim(j).twistAngle   )
end
subplot(325); xlabel('Time')
subplot(326); xlabel('Stroke Angle')
subplot(321); ylabel('In-plane angle','Rotation',0)
subplot(323); ylabel('Out-plane angle','Rotation',0)
subplot(325); ylabel('Twist angle','Rotation',0)
%% animations for check 

% % % figure()
% % % for j = 1:81
% % %    plotMesh( squeeze(sim(1).xyzPoints(j,:,:))' ,dotStyle)  
% % % %    plotMesh( squeeze(sim(2).xyzPoints(j,:,:))' ,dotStyle)  
% % %    axis([-5500 5500,-5500,5500,-5500,5500])
% % %    view(60,20)
% % %    drawnow
% % %    pause(0.1)
% % %    hold off
% % % end
% 
% %     end
% %
% % % figure()
% % % for j = 1:81
% % %    plotMesh( squeeze(sim(1).xyzHaltereFrame(j,:,:))' ,dotStyle)  
% % %    hold on 
% % %    plotMesh( squeeze(sim(2).xyzHaltereFrame(j,:,:))' ,dotStyle)  
% % %    axis([0,5000,-200,200,-300,300])
% % %    view(80,10)
% % %    hold off
% % %    drawnow
% % % %    pause(0.1)
% % % end
% %% 
% % % figure()
% % % for j = 1:81%n_times
% % %     scatter( squeeze(sim(1).xyzHaltereFrameMiddle(j,[2])), squeeze(sim(1).xyzHaltereFrameMiddle(j,[3])),'r'); hold on
% % %     scatter(  squeeze(sim(1).xyzHaltereFrame(j,sidePoints,2)), squeeze(sim(1).xyzHaltereFrame(j,sidePoints,3)),'k'); hold on
% % % 
% % %     scatter( squeeze(sim(2).xyzHaltereFrameMiddle(j,[2])), squeeze(sim(1).xyzHaltereFrameMiddle(j,[3])),'r'); hold on
% % %     scatter(  squeeze(sim(2).xyzHaltereFrame(j,sidePoints,2)), squeeze(sim(1).xyzHaltereFrame(j,sidePoints,3)),'k'); hold on
% % % 
% % %    axis([-200,200,-300,300])
% % % %    view(80,10)
% % %    drawnow
% % % %    pause(0.1)
% % %    hold off
% % % end
    

