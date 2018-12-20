clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
% loadName = 'figure2_deformMesh';
% saveName = 'figure2_deformMesh';

loadName = 'figure2_deformLines';
saveName = 'figure2_deformLines';
% renew_data_load = true
renew_data_load = false 
if renew_data_load
    FEA(1).name = 'Haltere_CraneFlyLowDensityWbulb_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_Sphere_Om10';
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
                eul_phi = euler_angle('Y',FEA(j).phi(k))^-1;
                eul_theta = euler_angle('Z',FEA(j).theta(k))^-1;
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
t_ind = 100;
deform_mult = 200;
OOP_mult = 2000;
xDes = [0:150:4650];

FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 ); 
FEA(1).xrtheta(:,3) = wrapTo2Pi(  atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi -0.02 )+0.02;


Xb3 = zeros(length(xDes), length(theta),200); 
Yb3 = Xb3;
Zb3 = Xb3; 
phiAngle = 0:0.005:1+1;
for j = 1:200
%     FEA(1).diffPerPoint = squeeze(FEA(1).xyzHaltereFrame(:,t_ind,:)) - FEA(1).xyz;
    j
%     if 1 ==2 
%      FEA(1).diffPerPoint(:,2) = - FEA(2).diffPerPoint(:,2);
%     end
    for k=1:length(xDes)
        dx = abs(FEA(1).xrtheta(:,1)-xDes(k) );
        dr = abs(FEA(1).xrtheta(:,2)-150 );
        for l = 1:length(theta)
    %         theta(k)
            da = abs(FEA(1).xrtheta(:,3) - theta(l) );
            J = dx.^2 + dr.^2+ (da*150.^2);
            [V,I] = min(J);
%             X1(k,l) = FEA(1).xyz(I,1)  + FEA(1).diffPerPoint(I,1)*deform_mult; 
%             Y1(k,l) = FEA(1).xyz(I,2) + FEA(1).diffPerPoint(I,2)*deform_mult; 
%             Z1(k,l) = FEA(1).xyz(I,3) + FEA(1).diffPerPoint(I,3)*deform_mult; 
%             C1(k,l) = FEA(1).strain(I,t_ind);
            
            Xb(k,l) = FEA(1).xyz(I,1) ; 
            Yb(k,l) = FEA(1).xyz(I,2) ; 
            Zb(k,l) = FEA(1).xyz(I,3) ; 
            
            
           Xb2(k,l,1) = FEA(1).xyzPoints(I,j,1);
           Yb2(k,l,1) = FEA(1).xyzPoints(I,j,2);
           Zb2(k,l,1) = FEA(1).xyzPoints(I,j,3);
           
           
            xyzTemp = [FEA(1).xyz(I,1), FEA(1).xyz(I,2), FEA(1).xyz(I,3) ];
%            for j2 = 1:length( 

% %             phiAngle(j) = pi/4;
%                 eul_1 = [ 1       0                     0;...
%                             0,  cos( phiAngle(j) ),  - sin(  phiAngle(j) )  ; ...
%                            0   sin(  phiAngle(j) ) cos(  phiAngle(j) )  ]^-1; 

            eul_1 = euler_angle('X',0)^-1;
            eul_2 = euler_angle('Y',phiAngle(j) )^-1;
            eul_3 = euler_angle('Z',0)^-1;
              xyzTranspose = eul_1*eul_2*eul_3*xyzTemp'; 
            Xb3(k,l,j) = xyzTranspose(1) ; 
            Yb3(k,l,j) = xyzTranspose(2) ; 
            Zb3(k,l,j) = xyzTranspose(3) ; 
            
           
        end
    end
end




%% Figure 2
fig2 = figure();
    width = 4;     % Width in inches,   find column width in paper 
    height = 4;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar

% subplot 311, rotation around x 
xc = 0; yc = 0; zc = 0; 
xr = 440; yr = 440; zr = 440; 
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

angles = [0,0.3,0];
eul_1 = euler_angle('X',angles(1))^-1;
eul_2 = euler_angle('Y',angles(2))^-1;
eul_3 = euler_angle('Z',angles(3))^-1;

for j = 1:size(x,1)
    for k = 1:size(x,2)
        xyzTemp = [x(j,k), y(j,k), z(j,k) ];
        xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
        xOm10(j,k) = xyzT(1);
        yOm10(j,k) = xyzT(2);
        zOm10(j,k) = xyzT(3);
    end
end
Cb = zeros(size(Xb));

Imax = 155; 
[z,y,x] = ellipsoid(0,0,0,440,440,440,16);

C = zeros(size(x));
hold on
% subplot(311); hold on 
% 	sb = surf(Xb,Yb,Zb,Cb);
%         set(sb,surfParamBackground{:})
% %     s2 = surf(Xj,Yj,Zj,Cj);
% %         set(s2,surfParamForeground{:})
%     s1 = surf(x+4800,y,z,C);
%         set(s1,surfParamBackground{:})
%         
%         
	s2 = surf(Xb2,Yb2,Zb2,Cb);
        set(s2,surfParamBackground{:})
        j = 200; 
	s3 = surf( Xb3(:,:,j),Yb3(:,:,j),Zb3(:,:,j),Cb);
        set(s3,surfParamBackground{:})
        
%     s3 = surf( xOm10 +4800+ FEA(2).diffPerPoint(Imax,1)*1.1*deform_mult ,...
%                 yOm10+ FEA(2).diffPerPoint(Imax,2)*1.1*deform_mult,...
%                 zOm10+ FEA(2).diffPerPoint(Imax,3)*1.1*deform_mult,...
%                 C);
%         set(s3,surfParamForeground{:}) 
        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z') 
        view(45,30)
        
        