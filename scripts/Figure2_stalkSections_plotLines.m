clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure2_deformLines';
saveName = 'figure2_deformLines';

renew_data_load = false
% renew_data_load = true
if renew_data_load 
    FEA(1).name = 'Haltere_CraneFlyLowDensityWbulb_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_Sphere_Om10'; 
    FEA(3).name = 'Haltere_CraneFlyLowDensitywBulb_sphereCrossStalk_Om0';
    FEA(4).name = 'Haltere_CraneFlyLowDensitywBulb_sphereCrossStalk_Om10';
    for j =  1:length(FEA)
        tic
        [~, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
        [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc 
    end
    
    er = 0.1;
    circleDistance = 5000;               % distance from base to haltere 
    bulbRadius = 500; 
    for j = 1:length(FEA)
        xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
        yMatch = find(   abs(abs( FEA(j).xyz(:,2) ) - bulbRadius) <er & ...
                     abs( FEA(j).xyz(:,3) )   <er);
        zMatch = find(   abs(abs( FEA(j).xyz(:,3) ) - bulbRadius)  <er & ...
                     abs( FEA(j).xyz(:,2) )   <er);
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
%% deformation in angles 

fig1 = figure();
    width = 3.5;     % Width in inches,   find column width in paper 
    height = 2.5;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
    lineType = {'-',':','--'};
    
for j = [1,3] 
    subplot(331); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) ,lineType{j} );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(334 ); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) ,lineType{j}   );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(337); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It)  ,lineType{j} );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end


%%  
for j = 2% 
    subplot(3,3,2 ); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(3,3,5 ); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(3,3,8 ); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end

%% Figure 2

axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15] ,...
               'YLim',[-1,1]*1e-4}; 
           
           
for j = 4 
    subplot(333); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(336); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(339); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end 
 

%% Setting paper size for saving 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
tightfig;
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_deformStalk' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_deformStalk'], '-dsvg', '-r600');
