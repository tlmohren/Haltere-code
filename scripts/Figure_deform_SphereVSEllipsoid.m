clc;clear all;close all
addpathFolderStructureHaltere()

loadName = 'figure3_data_FEA';
saveName = 'figure3_data_FEA_eval3000';
renew_data_load = true
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidHor_Om0';
    FEA(4).name = 'Haltere_CraneFly_ellipsoidHor_Om10';
    FEA(5).name = 'Haltere_CraneFly_ellipsoidVer_Om0';
    FEA(6).name = 'Haltere_CraneFly_ellipsoidVer_Om10';
    parameters = { 'disp','u2','v2','w2','flapangle','theta_angle', 'X (µm)'}; % select which parameters to load 
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, FEA(j).data, ~] = loadCSV( ['data' filesep  FEA(j).name], parameters);
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
        FEA(j).phi = FEA(j).data(1,:,5) ;
        FEA(j).theta = -FEA(j).data(1,:,6) ; 
        [ n_points,n_times, n_deform] = size( FEA(j).data(:,:,[2,3,4]) );
        % compute actual point locations 
        FEA(j).xyzPoints = FEA(j).data(:,:,[2,3,4])  + ...       
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
show_animation = false 

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

%%  deformation in micrometers

figure()
labels = {'$\Delta x$','$\Delta y$','$\Delta z$','$\gamma$'};
for j = 1:length(FEA)/2
    subplot(2,2,1); hold on 
        plot( FEA(j*2-1).twistAngle)
        ylabel( labels{4} );
    subplot(2,2,3); hold on 
        plot( FEA(j*2).twistAngle)
        ylabel( labels{4} );
    subplot(2,2,2); hold on
        sig =FEA(j*2-1).twistAngle(end-fftLen+1:end);
        [kf,y] = quick_fft(sig,1/dt*mult);
        plot(kf,y,'*-')
    subplot(2,2,4); hold on
        sig =FEA(j*2).twistAngle(end-fftLen+1:end);
        [kf,y] = quick_fft(sig,1/dt*mult);
        plot(kf,y,'*-')
end
subplot(224)
legend({'Sphere','Hor ellipse','Ver ellipse'})
