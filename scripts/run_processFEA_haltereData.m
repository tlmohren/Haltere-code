% process FEA results

loadName = 'FEA_processed_data';
saveName = 'FEA_processed_data';

renew_data_load = true
% renew_data_load = false

er = 0.1;
circleDistance = 5000;               % distance from base to haltere 
bulbA = [500,500, 500,500, 353.543,353.543,    1000,1000,     500,500,   353.543,353.543,   1000,1000,   ] ; 
bulbB = [500,500, 500,500, 1000,1000,          353.543,353.54,500,500,   1000,1000,    353.543,353.543,] ; 
yOffset = [0,0,0,0,0,0,0,0, 0,0, 0,0, 150, 150]; 
zOffset = [0,0,0,0,0,0,0,0,  150, 150, 150, 150, 0,0]; 
    
if renew_data_load
    
    % names of all simulations 
    FEA(1).name = 'Haltere_CraneFlyLowDensityWbulb_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_Sphere_Om10';
    FEA(3).name = 'Haltere_CraneFlyLowDensitywBulb_sphereCrossStalk_Om0';
    FEA(4).name = 'Haltere_CraneFlyLowDensitywBulb_sphereCrossStalk_Om10';
    FEA(5).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidVer_Om0';
    FEA(6).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_ellipsoidVer_Om10';
    FEA(7).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidHor_Om0';
    FEA(8).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_ellipsoidHor_Om10';
    FEA(9).name =  'Haltere_CraneFlyLowDensitywBulb_SphereVerOffset_Om0';
    FEA(10).name =  'Haltere_CraneFlyLowDensitywBulb_SphereVerOffset_Om10';
    FEA(11).name =  'Haltere_CraneFlyLowDensitywBulb_ellipsoidVerOffset_Om0';
    FEA(12).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidVerOffset_Om10';   
    FEA(13).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidHorOffset_Om0';
    FEA(14).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidHorOffset_Om10';

   % load different properties from the data 
    for j =  1:length(FEA)
        tic
        [~, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
        [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc 
    end
    
    % find points for top and side of haltere, (need only size? 
    for j = 1:length(FEA)
        xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
        yMatch = find(   abs(abs( FEA(j).xyz(:,2) -yOffset(j)) - bulbA(j)   ) <er & ...
                     abs( abs( FEA(j).xyz(:,3)-zOffset(j) )   ) <er) ;
        zMatch = find(   abs(abs( FEA(j).xyz(:,3)-zOffset(j) ) - bulbB(j)   ) <er & ...
                     abs( abs( FEA(j).xyz(:,2)-yOffset(j) )   )<er );
        FEA(j).sideInds = intersect(xMatch,yMatch);
%         FEA(j).topInds = intersect(xMatch,zMatch);
    end
    
    % Transpose deformation to haltere frame 
    for j =  1:length(FEA)
        tic
        % rename angles 
        FEA(j).phi = FEA(j).angles(1,:,1) ;
        FEA(j).theta = -FEA(j).angles(1,:,2) ; 
        [ n_points,n_times, n_deform] = size( FEA(j).deform ); 
        
        % compute actual point locations 
        FEA(j).xyzPoints = FEA(j).deform  + ...       
            permute( repmat( FEA(j).xyz,1,1, n_times), [1,3,2] ) ;
        
        % transpose moving/deforming mesh points into rigid haltere frame 
        for k = 1:n_times
            for l = 1:n_points
                eul_phi = euler_angle('Y',FEA(j).phi(k))^-1;
                eul_theta = euler_angle('Z',FEA(j).theta(k))^-1;
                FEA(j).xyzHaltereFrame(l,k,:) = eul_phi* eul_theta*  squeeze(FEA(j).xyzPoints(l,k,:) ) ;
            end
        end
    %     compute middle of haltere bulb 
        FEA(j).xyzHaltereFrameMiddle = squeeze(  mean( FEA(j).xyzHaltereFrame( FEA(j).sideInds,:,:), 1)  )- ones(length(FEA(j).phi),1)*[0,yOffset(j),zOffset(j)];

        % compute actual difference between two 
        FEA(j).dz = diff( FEA(j).xyzHaltereFrame( FEA(j).sideInds,:,3));
        FEA(j).dy = diff( FEA(j).xyzHaltereFrame( FEA(j).sideInds,:,2));
        
        % compute three deformation angles 
        FEA(j).twistAngle = atan2(FEA(j).dz, abs( FEA(j).dy) );
        FEA(j).yAngle = atan( FEA(j).xyzHaltereFrameMiddle(:,3) ./ FEA(j).xyzHaltereFrameMiddle(:,1));
        FEA(j).zAngle = atan( FEA(j).xyzHaltereFrameMiddle(:,2) ./ FEA(j).xyzHaltereFrameMiddle(:,1));
        toc 
    end
    
    for j = 1:2
        circleDistance = 300;               % distance from base to haltere 
        circleRadius = 150;                 % radius of haltere   
        mindist =  min( abs( FEA(j).xyz(:,1) - circleDistance) );
        xMatch = find(  abs(FEA(j).xyz(:,1) - circleDistance) <= (mindist+1) );
        rMatch = find( sqrt(FEA(j).xyz(:,2).^2 + FEA(j).xyz(:,3).^2)  >= circleRadius*0.99 );

        yMatch = find( round( abs( FEA(j).xyz(:,2) ), 7) == circleRadius );
        zMatch = find( round( abs( FEA(j).xyz(:,3) ), 7) == circleRadius );

        FEA(j).circleIndsUnsorted = intersect( xMatch,rMatch); 
        
        angle = atan2( FEA(j).xyz( FEA(j).circleIndsUnsorted,3), ...
            FEA(j).xyz( FEA(j).circleIndsUnsorted,2) );
        angleDeg = rad2deg(angle);
        angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
        [V,I_sort] = sort(angleDeg,'ascend');
% 
        FEA(j).circleInds= FEA(j).circleIndsUnsorted(I_sort);
    end
    
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end

FEA