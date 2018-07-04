clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure5_deform';
saveName = 'figure5_deform';

renew_data_load = false
% renew_data_load = true
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om10';
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
        [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc 
    end
    
    sideL = 243.435; 
    circleDistance = 3000;               % distance from base to haltere 
    er = 0.01;
    for j = 1:length(FEA)
        xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
        yMatch = find(   abs(abs( FEA(j).xyz(:,2) ) - sideL) <er & ...
                     abs( FEA(j).xyz(:,3) )   <er);
        zMatch = find(   abs(abs( FEA(j).xyz(:,3) ) - sideL)  <er & ...
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
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

deformLabels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};

axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-0.02,0.02]}; 
axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-1,1]*1e-3}; 
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*1e-3}; 

for j = [1]
    subplot(3,1,1); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(312); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
%         
    subplot(313); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end

%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure5_deformAnglePlotOM0' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure5_deformAnglePlotOM0'], '-dsvg', '-r600');


%%



fig2 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size

len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

deformLabels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};

axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-0.02,0.02]}; 
axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-1,1]*1e-3}; 
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*1e-3}; 

for j = 2
    subplot(3,1,1); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(312); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
%         
    subplot(313); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end

%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure5_deformAnglePlotOM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure5_deformAnglePlotOM10'], '-dsvg', '-r600');