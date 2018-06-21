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


%% 
show_animation = false;

if show_animation 
    dotStyle = {'Marker','.','MarkerFaceColor','k'};
    dotStyle2 = {'Marker','.','MarkerFaceColor','r'};
    figure('Position',[300,300,800,400])
    for j = 1:100
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

%% deformation in angles 
dt = 0.001;
mult = 1;
fftLen = 200;
t = 0:0.001:0.35;
labels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};

figure(); hold on 
for j = 1:length(FEA)/2
    subplot(3,2,1); hold on 
        plot(t,FEA(j*2-1).yAngle )
        plot(t,FEA(j*2).yAngle )
        xlabel('Time (s)'); ylabel( labels{1} );
    subplot(3,2,2); hold on 
        plot(FEA(j).phi, FEA(j*2-1).yAngle )
        plot(FEA(j).phi, FEA(j*2).yAngle  )
        xlabel('$\phi$'); ylabel( labels{1} );
    subplot(323); hold on 
        plot(t,FEA(j*2-1).zAngle )
        plot(t,FEA(j*2).zAngle )
        xlabel('Time (s)'); ylabel( labels{2} );
    subplot(324); hold on 
        plot(FEA(j).phi, FEA(j*2-1).zAngle )
        plot(FEA(j).phi, FEA(j*2).zAngle  )
        xlabel('$\phi$'); ylabel( labels{2} );
    subplot(325); hold on 
        plot(t,FEA(j*2-1).twistAngle)
        plot(t,FEA(j*2).twistAngle )
        xlabel('Time (s)'); ylabel( labels{3} );
    subplot(326); hold on 
        plot(FEA(j).phi, FEA(j*2-1).twistAngle )
        plot(FEA(j).phi, FEA(j*2).twistAngle)
        xlabel('$\phi$'); ylabel( labels{3} );
end

%% deformation in micrometers 
labels = {'$\Delta x$','$\Delta y$','$\Delta z$','$\gamma$'};
for j = 1:length(FEA)/2
figure()
    for k = 1:3
    subplot(4,2,k*2-1); hold on 
        plot(FEA(j*2-1).xyzHaltereFrameMiddle(:,k))
        plot(FEA(j*2).xyzHaltereFrameMiddle(:,k))
        ylabel( labels{k} );
    
    subplot(4,2,k*2); hold on
        sig = FEA(j*2-1).xyzHaltereFrameMiddle(end-fftLen+1:end,k);
        [kf,y] = quick_fft(sig,1/dt*mult);
        plot(kf,y,'*-')
        sig = FEA(j*2).xyzHaltereFrameMiddle(end-fftLen+1:end,k);
        [kf,y] = quick_fft(sig,1/dt*mult);
        plot(kf,y,'*-')
    end
    subplot(427);hold on 
        plot( FEA(j*2-1).twistAngle)
        plot( FEA(j*2).twistAngle)
        ylabel( labels{4} );
    subplot(428); hold on
        sig =FEA(j*2-1).twistAngle(end-fftLen+1:end);
        [kf,y] = quick_fft(sig,1/dt*mult);
        plot(kf,y,'*-')
        sig =FEA(j*2).twistAngle(end-fftLen+1:end);
        [kf,y] = quick_fft(sig,1/dt*mult);
        plot(kf,y,'*-')
end

%% 

    %
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
OOP_mult = 300;
xDes = [0:150:4800];
surfParamBackground = {'FaceAlpha',0.2,'EdgeAlpha',0.2};
surfParamForeground = {'EdgeAlpha',0.2};
    
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

xc = 0;
zc = 0;
yc = 0; 
xr = 300; 
yr = 300;
zr = 946; 
n = 16;

[z,y,x] = ellipsoid(zc,yc,xc,zr,yr,xr,n);

angles = [ax, ay, az];

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
Cb = zeros(size(Xb));

[zb,yb,xb] = ellipsoid(0,0,0,948,300,300,16);

figure('Position',[100 100 1200 500]); hold on
subplot(121); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackground{:})
        c_ell = zeros(size(xb));
    s1 = surf(xb+5e3,yb,zb,c_ell);
        set(s1,surfParamBackground{:})
    s2 = surf(Xj,Yj,Zj,Cj);
        set(s2,surfParamForeground{:})
        colormap(strainScheme)%     colorbar
        C = zeros(size(x));
        set(s1,surfParamBackground{:})
        Imax = 155;
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

angles = [dax, day, daz];
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
        xOMdiff(j,k) = xyzT(1);
        yOMdiff(j,k) = xyzT(2);
        zOMdiff(j,k) = xyzT(3);
    end
end

        [zb,yb,xb] = ellipsoid(0,0,0,948,300,300,16);
subplot(122); hold on 
	sb = surf(Xb,Yb,Zb,Cb);
        set(sb,surfParamBackground{:})
        c_ell = zeros(size(xb));
    s1 = surf(xb+5e3,yb,zb,c_ell);
        set(s1,surfParamBackground{:})
    s2 = surf(XjDiff,YjDiff,ZjDiff,CjDiff);
        set(s2,surfParamForeground{:})
        colormap(strainScheme)%     colorbar
        C = zeros(size(x));
        Imax = 155;
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
