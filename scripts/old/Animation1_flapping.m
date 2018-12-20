clc;clear all;close all
addpathFolderStructureHaltere()
%  
loadName = 'figure4_strainData';
saveName = 'figure4_strainData';
renew_data_load = false
% renew_data_load = true

if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    for j =  1:length(FEA)
        tic
        [~, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
        [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc 
    end
    % Determine Circle locations
    for j = 1:length(FEA)
        circleDistance = 3000;               % distance from base to haltere 
        circleRadius = 150;                 % radius of haltere   
        mindist =  min( abs( FEA(j).xyz(:,1) - circleDistance) );
        xMatch = find(  abs(FEA(j).xyz(:,1) - circleDistance) <= (mindist+1) );
        yMatch = find( round( abs( FEA(j).xyz(:,2) ), 7) == circleRadius );
        zMatch = find( round( abs( FEA(j).xyz(:,3) ), 7) == circleRadius );
        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    % Transpose deformation to haltere frame 
    for j = 1:length(FEA)
        tic
        FEA(j).phi = FEA(j).angles(1,:,1) ;
        FEA(j).theta = -FEA(j).angles(1,:,2) ; 
        [ n_points,n_times, n_deform] = size( FEA(j).deform );
        % compute actual point locations 
        FEA(j).xyzPoints = FEA(j).deform  + ...       
            permute( repmat( FEA(j).xyz,1,1, n_times), [1,3,2] ) ;
        % transpose points into rigid haltere frame 
        for k = 1:n_times
            for l = 1:n_points
                eul_phi = [cos( FEA(j).phi(k))      0       sin( FEA(j).phi(k)) ; ...
                            0              1       0 ;...
                            -sin( FEA(j).phi(k))    0        cos( FEA(j).phi(k))]^-1;

                eul_theta = [cos( FEA(j).theta(k))   sin( FEA(j).theta(k))    0 ; ...
                             -sin( FEA(j).theta(k)) cos( FEA(j).theta(k))     0 ;...
                             0          0               1]^-1;

                FEA(j).xyzHaltereFrame(l,k,:) = eul_phi* eul_theta*  squeeze(FEA(j).xyzPoints(l,k,:) ) ;
            end
        end
    %     compute middle point j
        FEA(j).xyzHaltereFrameMiddle = squeeze(  mean( FEA(j).xyzHaltereFrame( FEA(j).sideInds,:,:), 1)  );

        FEA(j).dz = diff( FEA(j).xyzHaltereFrame( FEA(j).sideInds,:,3));
        FEA(j).dy = diff( FEA(j).xyzHaltereFrame( FEA(j).sideInds,:,2));
        FEA(j).twistAngle = atan2(FEA(j).dz, abs( FEA(j).dy) );
        FEA(j).yAngle = atan( FEA(j).xyzHaltereFrameMiddle(:,3) ./ FEA(j).xyzHaltereFrameMiddle(:,1));
        FEA(j).zAngle = atan( FEA(j).xyzHaltereFrameMiddle(:,2) ./ FEA(j).xyzHaltereFrameMiddle(:,1));
        toc 
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end
%% 
if 0
    dotStyle = {'Marker','.','MarkerFaceColor','k'};
    dotStyle2 = {'Marker','.','MarkerFaceColor','r'};
    figure('Position',[300,300,800,400])
    for j = 90:100
        subplot(121)
           title( ['t = ' num2str(j/1e3)])
           plotMesh( squeeze(FEA(1).xyzPoints(:,j,:))' ,dotStyle)  
           plotMesh( squeeze(FEA(2).xyzPoints(:,j,:))' ,dotStyle2)  
           axis([-5500 5500,-5500,5500,-5500,5500])
           view(60,20)
           drawnow
           pause(0.01)
           hold off
       subplot(122)
       plotMesh( squeeze(FEA(1).xyzHaltereFrame(:,j,:))' ,dotStyle)  
       plotMesh( squeeze(FEA(2).xyzHaltereFrame(:,j,:))' ,dotStyle2)  
           axis([0,5000,-200,200,-300,300])
           view(80,10)
           hold off
           drawnow
    end
end

%% 
redPurple = [158,1,66
213,62,79
244,109,67
253,174,97
254,224,139
255,255,191
230,245,152
171,221,164
102,194,165
50,136,189
94,79,16];
strainScheme = colorSchemeInterp(redPurple/255, 500);

%%
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;

t_ind = 100;
deform_mult = 30;
OOP_mult = 300;
xDes = [0:150:4800];
% surfParamBackground = {'FaceAlpha',0.2,'EdgeAlpha',0.2};
% surfParamForeground = {'EdgeAlpha',0.2};
    
surfParamBackground = {'FaceAlpha',0.2,'EdgeAlpha',0.1};
surfParamBackgroundTwist = {'FaceAlpha',0.2,'EdgeAlpha',0.2};
surfParamBackgroundTwistStalk = {'FaceAlpha',0,'EdgeAlpha',0.1};
surfParamForeground = {'EdgeAlpha',0.4};
%% 
FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
FEA(1).xrtheta(:,3) = wrapTo2Pi(  atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi -0.02 )+0.02;
FEA(2).xrtheta = FEA(1).xrtheta;


for j = 1:200
    for k=1:length(xDes)
        dx = abs(FEA(1).xrtheta(:,1)-xDes(k) );
        dr = abs(FEA(1).xrtheta(:,2)-150 );
        for l = 1:length(theta)
    %         theta(k)
            da = abs(FEA(1).xrtheta(:,3) - theta(l) );
            J = dx.^2 + dr.^2+ (da*150.^2);
            [V,I] = min(J);
           Xb(k,l,j) = FEA(1).xyzPoints(I,j,1);
           Yb(k,l,j) = FEA(1).xyzPoints(I,j,2);
           Zb(k,l,j) = FEA(1).xyzPoints(I,j,3);
           Cb(k,l,j) = FEA(1).strain(I,j);
%            k,l
        end
    end
end
FEA_n = 1; 

ax = FEA(FEA_n ).twistAngle(t_ind)*deform_mult*1.5;
ay = FEA(FEA_n ).yAngle(t_ind)*deform_mult*1.5;
az = FEA(FEA_n ).zAngle(t_ind)*deform_mult*1.5;

xc = 0;
zc = 0;
yc = 0; 
xr = 440; 
yr = 440;
zr = 440; 
n = 16;

[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

[zb,yb,xb] = ellipsoid(0,0,0,440,440,440,16);

%0-----------------------------------
fig2 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar


for frame = 100:200
    plot3( [0,1]*4e3,[0,1]*0,[0,1]*0 ,'k') 
    hold on
    plot3( [0,1]*0,[0,-1]*4e3,[0,1]*0 ,'k') 
    plot3( [0,1]*0,[0,1]*0,[0,1]*4e3 ,'k') 
     text(4.3e3,0,0,'y')
     text(0,-5.3e3,0,'x')
     text(-100,0,4.5e3,'z')
    angles = [ax, -FEA(FEA_n ).phi(frame), az];

    for j = 1:size(x,1)
        for k = 1:size(x,2)
            xyzTemp = [x(j,k), y(j,k), z(j,k) ];

                    eul_1 = [ 1       0                     0;...
                                0,  cos( angles(1) ),  - sin(  angles(1) )  ; ...
                               0   sin(  angles(1) ) cos(  angles(1) )  ]^-1;
                    eul_2 = [cos(  angles(2))      0       sin( angles(2)) ; ...
                                0              1       0 ;...
                                -sin( angles(2))    0        cos(  angles(2))]^-1;

                    eul_3 = [cos(  angles(3))   sin(  angles(3))    0 ; ...
                                 -sin(  angles(3)) cos(  angles(3))     0 ;...
                                 0          0               1]^-1;
            xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
            xOm10(j,k) = xyzT(1);
            yOm10(j,k) = xyzT(2);
            zOm10(j,k) = xyzT(3);
        end
    end


    s2 = surf(Xb(:,:,frame),Yb(:,:,frame),Zb(:,:,frame),Cb(:,:,frame));
        set(s2,surfParamForeground{:})
        
        colormap(strainScheme)%     colorbar
        C = zeros(size(x));
        ImaxBottom = 155;
        ImaxTop = 39;
        centerPoint = [( FEA(FEA_n ).xyzPoints(ImaxBottom,frame,1) + FEA(FEA_n ).xyzPoints(ImaxTop,frame,1))/2;
                ( FEA(FEA_n ).xyzPoints(ImaxBottom,frame,2) + FEA(FEA_n ).xyzPoints(ImaxTop,frame,2))/2;
                ( FEA(FEA_n ).xyzPoints(ImaxBottom,frame,3) + FEA(FEA_n ).xyzPoints(ImaxTop,frame,3))/2];

    % bulb angle 
    s3 = surf( xOm10 + centerPoint(1),...
                yOm10 + centerPoint(2) ,...
                zOm10 + centerPoint(3)  ,...
                C);
        set(s3,surfParamForeground{:})
        axis equal
        caxis([-1,1]*0.0005)
        axis([-300 5500 -5500 500 -5500 5500])
        shading faceted
        view(30,15)
        axis off
        drawnow
        hold off 

end