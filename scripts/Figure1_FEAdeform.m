clc;clear all;close all
addpathFolderStructureHaltere()

loadName = 'figureS2_deformMesh';
saveName = 'figureS2_deformMesh';
renew_data_load = false
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
%     parameters = { 'disp','u2','v2','w2','flapangle','theta_angle', 'X (µm)'}; % select which parameters to load 
    for j =  1:length(FEA)
        tic
%         [FEA(j).xyz, FEA(j).data, ~] = loadCSV( ['data' filesep  FEA(j).name], parameters);
        [~, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
%         [FEA(j).xyz, FEA(j).a, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
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

%% deformation in angles 

width = 2;     % Width in inches,   find column width in paper 
height = 3;    % Height in inches
fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

t = 0:0.001:0.35;
labels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};
axOpts1 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 
axOpts2 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 
axOpts3 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 

for j = 1:length(FEA)/2
    subplot(3,1,1); hold on 
        plot(t,FEA(j*2-1).yAngle )
        plot(t,FEA(j*2).yAngle )
        ylabel( labels{1} );
        ax = gca();
        set(ax,axOpts1{:})
    subplot(312); hold on 
        plot(t,FEA(j*2-1).zAngle )
        plot(t,FEA(j*2).zAngle )
        ylabel( labels{2} );
        ax = gca();
        set(ax,axOpts2{:})
    subplot(313); hold on 
        plot(t,FEA(j*2-1).twistAngle)
        plot(t,FEA(j*2).twistAngle )
        xlabel('Time (s)'); ylabel( labels{3} );
        ax = gca();
        set(ax,axOpts3{:})
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
% 
strainScheme = colorSchemeInterp(redPurple/255, 500);


%%
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;

t_ind = 100;
deform_mult = 30;
OOP_mult = -300;
xDes = [0:150:4800];
surfParamBackground = {'FaceAlpha',0.2,'EdgeAlpha',0.1};
surfParamBackgroundTwist = {'FaceAlpha',0.4,'EdgeAlpha',0.3};
surfParamBackgroundTwistStalk = {'FaceAlpha',0,'EdgeAlpha',0.1};
surfParamForeground = {'EdgeAlpha',0.4};
    
%% 
FEA(1).xrtheta(:,1) = FEA(1).xyz(:,1);
FEA(1).xrtheta(:,2) = sqrt( FEA(1).xyz(:,2).^2  +  FEA(1).xyz(:,3).^2 );
FEA(1).xrtheta(:,3) = atan2( FEA(1).xyz(:,3), FEA(1).xyz(:,2) ) +pi ;
FEA(2).xrtheta = FEA(1).xrtheta;

for j = 1:length(FEA)
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


        
        

%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_deformAnglePlot' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_deformAnglePlot'], '-dsvg', '-r600');
























%% Figure 2
width = 2;     % Width in inches,   find column width in paper 
height = 3;    % Height in inches
fig2 = figure();
set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
colormap(strainScheme)%     colorbar

%% subplot 311, rotation around x 
xc = 0;
zc = 0;
yc = 0; 
xr = 300; 
yr = 300;
zr = 946; 
n = 16;
[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

angles = [0,0.1,0];
eul_1 = [ 1       0                     0;...
            0,  cos( angles(1) ),  - sin(  angles(1) )  ; ...
           0   sin(  angles(1) ) cos(  angles(1) )  ]^-1;
eul_2 = [cos(  angles(2))      0       sin( angles(2)) ; ...
            0              1       0 ;...
            -sin( angles(2))    0        cos(  angles(2))]^-1;
eul_3 = [cos(  angles(3))   sin(  angles(3))    0 ; ...
             -sin(  angles(3)) cos(  angles(3))     0 ;...
             0          0               1]^-1;
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

[z,y,x] = ellipsoid(0,0,0,948,300,300,16);
subplot(311); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackground{:})
        c_ell = zeros(size(x));
    s2 = surf(Xj,Yj,Zj,Cj);
        set(s2,surfParamForeground{:})
        C = zeros(size(x));
        Imax = 155;
    s1 = surf(x+5e3,y,z,c_ell);
        set(s1,surfParamBackground{:})
        
    s3 = surf( xOm10 +5000+ FEA(2).diffPerPoint(Imax,1)*deform_mult ,...
                yOm10+ FEA(2).diffPerPoint(Imax,2)*deform_mult,...
                zOm10+ FEA(2).diffPerPoint(Imax,3)*deform_mult,...
                C);
        set(s3,surfParamForeground{:})
        axis off
        axis equal
        axis tight;
        shading faceted
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        %% 
angles = [0,0,0.1];
eul_1 = [ 1       0                     0;...
            0,  cos( angles(1) ),  - sin(  angles(1) )  ; ...
           0   sin(  angles(1) ) cos(  angles(1) )  ]^-1;
eul_2 = [cos(  angles(2))      0       sin( angles(2)) ; ...
            0              1       0 ;...
            -sin( angles(2))    0        cos(  angles(2))]^-1;
eul_3 = [cos(  angles(3))   sin(  angles(3))    0 ; ...
             -sin(  angles(3)) cos(  angles(3))     0 ;...
             0          0               1]^-1;
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
        colormap(strainScheme)%     colorbar
    s1 = surf(x+5e3,y,z,c_ell);
        set(s1,surfParamBackground{:})
    s3 = surf( xOMdiff +5000+ totDiff(Imax ,1)*1.1 *OOP_mult ,...
            yOMdiff+ totDiff(Imax ,2)*1.1*OOP_mult,...
            zOMdiff+ totDiff(Imax ,3)*1.1*OOP_mult,...
            C);
        set(s3,surfParamForeground{:})

        axis tight;
        shading faceted
        axis off
        axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        
        
        
%% twist 
angles = [0.8,0,0];
eul_1 = [ 1       0                     0;...
            0,  cos( angles(1) ),  - sin(  angles(1) )  ; ...
           0   sin(  angles(1) ) cos(  angles(1) )  ]^-1;
eul_2 = [cos(  angles(2))      0       sin( angles(2)) ; ...
            0              1       0 ;...
            -sin( angles(2))    0        cos(  angles(2))]^-1;
eul_3 = [cos(  angles(3))   sin(  angles(3))    0 ; ...
             -sin(  angles(3)) cos(  angles(3))     0 ;...
             0          0               1]^-1;
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
    eul_1 = [ 1       0                     0;...
                0,  cos( angleTemp(1) ),  - sin(  angleTemp(1) )  ; ...
               0   sin(  angleTemp(1) ) cos(  angleTemp(1) )  ]^-1;
    eul_2 = [cos(  angleTemp(2))      0       sin( angleTemp(2)) ; ...
                0              1       0 ;...
                -sin( angleTemp(2))    0        cos(  angleTemp(2))]^-1;
    eul_3 = [cos(  angleTemp(3))   sin(  angleTemp(3))    0 ; ...
                 -sin(  angleTemp(3)) cos(  angleTemp(3))     0 ;...
                 0          0               1]^-1;
         
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
    s2 = surf(Xtwist,Ytwist,Ztwist,-Ctwist);
        set(s2,surfParamForeground{:})
        colormap(strainScheme)%     colorbar
    s1 = surf(x+5e3,y,z,c_ell);
        set(s1,surfParamBackgroundTwist{:})
    s3 = surf( xOMdiff +5000 ,...
            yOMdiff,...
            zOMdiff,...
            C);
        set(s3,surfParamForeground{:})

        axis tight;
        shading faceted
        axis off
        axis equal
        xlabel('x');ylabel('y');zlabel('z')
        view(40,40)
        

        
        
        
        
        
        
        
        

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
        
