clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure1_deform';
saveName = 'figure1_deform';

renew_data_load = false
% renew_data_load = true
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
    FEA(4).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
    FEA(5).name = 'Haltere_CraneFly_ellipsoidVer_Om0';
    FEA(6).name = 'Haltere_CraneFly_ellipsoidVer_Om10';   
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
deform_mult = 40;
OOP_mult = -600;
xDes = [0:150:4800];

FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
% FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;
FEA(1).xrtheta(:,3) = wrapTo2Pi(  atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi -0.02 )+0.02;


FEA(2).xrtheta = FEA(1).xrtheta;

for j = 1:2%length(FEA)
%     FEA(j).xrtheta = FEA(1).xrtheta;
    FEA(j).diffPerPoint = squeeze(FEA(j).xyzHaltereFrame(:,t_ind,:)) - FEA(j).xyz;
    for k=1:length(xDes)
        dx = abs(FEA(j).xrtheta(:,1)-xDes(k) );
        dr = abs(FEA(j).xrtheta(:,2)-150 );
        for l = 1:length(theta)
    %         theta(k)
            da = abs(FEA(j).xrtheta(:,3) - theta(l) );
            J = dx.^2 + dr.^2+ (da*150.^2);
            [V,I] = min(J);
            Xj(k,l) = FEA(j).xyz(I,1)  + FEA(j).diffPerPoint(I,1)*deform_mult; 
            Yj(k,l) = FEA(j).xyz(I,2) + FEA(j).diffPerPoint(I,2)*deform_mult; 
            Zj(k,l) = FEA(j).xyz(I,3) + FEA(j).diffPerPoint(I,3)*deform_mult; 
            Cj(k,l) = FEA(j).strain(I,t_ind);
            
            Xb(k,l) = FEA(j).xyz(I,1) ; 
            Yb(k,l) = FEA(j).xyz(I,2) ; 
            Zb(k,l) = FEA(j).xyz(I,3) ; 
        end
    end
end

totDiff = [(FEA(2).diffPerPoint(:,1)-FEA(1).diffPerPoint(:,1)),...
            (FEA(2).diffPerPoint(:,2)-FEA(1).diffPerPoint(:,2)), ...
            (FEA(2).diffPerPoint(:,3)-FEA(1).diffPerPoint(:,3))];

for k=1:length(xDes)
    dx = abs(FEA(j).xrtheta(:,1)-xDes(k) );
    dr = abs(FEA(j).xrtheta(:,2)-150 );
    for l = 1:length(theta)
        da = abs(FEA(j).xrtheta(:,3) - theta(l) );
        J = dx.^2 + dr.^2+ (da*150.^2);
        [V,I] = min(J);
        
        XjDiff(k,l) = FEA(2).xyz(I,1)  + (FEA(2).diffPerPoint(I,1)-FEA(1).diffPerPoint(I,1)) *OOP_mult ; 
        YjDiff(k,l) = FEA(2).xyz(I,2) + (FEA(2).diffPerPoint(I,2)-FEA(1).diffPerPoint(I,2)) *OOP_mult ; 
        ZjDiff(k,l) = FEA(2).xyz(I,3)+ (FEA(2).diffPerPoint(I,3)-FEA(1).diffPerPoint(I,3)) *OOP_mult ; 
        CjDiff(k,l) = FEA(2).strain(I,t_ind)-FEA(1).strain(I,t_ind);
    end
end

FEA(2).diffPerPointDiff = squeeze(FEA(j).xyzHaltereFrame(:,t_ind,:)) - FEA(j).xyz;
ax = FEA(2).twistAngle(t_ind)*deform_mult*1.5;
ay = FEA(2).yAngle(t_ind)*deform_mult*1.5;
az = FEA(2).zAngle(t_ind)*deform_mult*1.5;

dax = (FEA(2).twistAngle(t_ind)-FEA(1).twistAngle(t_ind)) *OOP_mult*1.85;
day = (FEA(2).yAngle(t_ind)-FEA(1).yAngle(t_ind)) *OOP_mult*1.85;
daz = (FEA(2).zAngle(t_ind)-FEA(1).zAngle(t_ind)) *OOP_mult*1.85;

surfParamBackground = {'FaceAlpha',0.2,'EdgeAlpha',0.1};
surfParamBackgroundTwist = {'FaceAlpha',0.2,'EdgeAlpha',0.2};
surfParamBackgroundTwistStalk = {'FaceAlpha',0,'EdgeAlpha',0.1};
surfParamForeground = {'EdgeAlpha',0.4};
    

%% Figure 2
fig2 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar

% subplot 311, rotation around x 
xc = 0; yc = 0; zc = 0;
xr = 300; yr = 300;zr = 946; 
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
[z,y,x] = ellipsoid(0,0,0,948,300,300,16);

C = zeros(size(x));

subplot(311); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackground{:})
    s2 = surf(Xj,Yj,Zj,Cj);
        set(s2,surfParamForeground{:})
    s1 = surf(x+5e3,y,z,C);
        set(s1,surfParamBackground{:})
    s3 = surf( xOm10 +5000+ FEA(2).diffPerPoint(Imax,1)*deform_mult ,...
                yOm10+ FEA(2).diffPerPoint(Imax,2)*deform_mult,...
                zOm10+ FEA(2).diffPerPoint(Imax,3)*deform_mult,...
                C);
        set(s3,surfParamForeground{:})

        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,17)
        
angles = [0,0,0.1];
eul_1 = euler_angle('X',angles(1))^-1;
eul_2 = euler_angle('Y',angles(2))^-1;
eul_3 = euler_angle('Z',angles(3))^-1;

for j = 1:size(x,1)
    for k = 1:size(x,2)
        xyzTemp = [x(j,k), y(j,k), z(j,k) ];
        xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
        xOMdiff(j,k) = xyzT(1);
        yOMdiff(j,k) = xyzT(2);
        zOMdiff(j,k) = xyzT(3);
    end
end


[z,y,x] = ellipsoid(0,0,0,948,300,300,16);
subplot(312); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackground{:})
    s2 = surf(XjDiff,YjDiff,ZjDiff,CjDiff);
        set(s2,surfParamForeground{:})
    s1 = surf(x+5e3,y,z,C);
        set(s1,surfParamBackground{:})
    s3 = surf( xOMdiff +5000+ totDiff(Imax ,1)*1.1 *OOP_mult ,...
            yOMdiff+ totDiff(Imax ,2)*1.1*OOP_mult,...
            zOMdiff+ totDiff(Imax ,3)*1.1*OOP_mult,...
            C);
        set(s3,surfParamForeground{:})
        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,17)
        
%twist 
angles = [0.8,0,0];
eul_1 = euler_angle('X',angles(1))^-1;
eul_2 = euler_angle('Y',angles(2))^-1;
eul_3 = euler_angle('Z',angles(3))^-1;

for j = 1:size(x,1)
    for k = 1:size(x,2)
        xyzTemp = [x(j,k), y(j,k), z(j,k) ];
        xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
        xOMdiff(j,k) = xyzT(1);
        yOMdiff(j,k) = xyzT(2);
        zOMdiff(j,k) = xyzT(3);
    end
end

for j = 1:size(Xb,1)
    angleTemp = angles*j/size(Xb,1);
    eul_1 = euler_angle('X',angleTemp(1))^-1;
    eul_2 = euler_angle('Y',angleTemp(2))^-1;
    eul_3 = euler_angle('Z',angleTemp(3))^-1;

    for k = 1:size(Xb,2)
        xyzTemp = [Xb(j,k), Yb(j,k), Zb(j,k) ];
        xyzT = eul_1*eul_2*eul_3*xyzTemp'; 
        Xtwist(j,k) = xyzT(1);
        Ytwist(j,k) = xyzT(2);
        Ztwist(j,k) = xyzT(3);
    end
end

Ctwist = ones(size(Cb));
Ctwist(1,1) = -1.5;
Ctwist(1,2) = 1.5;
[z,y,x] = ellipsoid(0,0,0,948,300,300,16);

subplot(313); hold on   
%     plot stalks 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackgroundTwistStalk{:})
    s2 = surf(Xtwist,Ytwist,Ztwist,Cb);
        set(s2,surfParamForeground{:})
        colormap(strainScheme)%     colorbar
    s1 = surf(x+5e3,y,z,C);
        set(s1,surfParamBackgroundTwist{:})
    s3 = surf( xOMdiff +5000 ,...
            yOMdiff,...
            zOMdiff,...
            C);
        set(s3,surfParamForeground{:})

        shading faceted
        axis tight;  axis off; axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,17)
        
%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
tightfig;
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure1_deformMesh' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure1_deformMesh'], '-dsvg', '-r600');
        
