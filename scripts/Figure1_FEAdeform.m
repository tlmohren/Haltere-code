clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure1_deform';
saveName = 'figure1_deform';

renew_data_load = false
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
%% deformation in angles 

fig1 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

deformLabels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};
        
axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15],'XTickLabel',{'','','',''} ,...
               'YLim',[-0.02,0.02]}; 
axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15],'XTickLabel',{'','','',''} ,...
               'YLim',[-1,1]*1e-3}; 
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15] ,...
               'YLim',[-1,1]*5e-5}; 
           
legend_entries = {'sphere','ellipsoid hor','ellipsoid ver'};

for j = [1,3,5]
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
        
    subplot(313); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end
    subplot(311)
        legend(legend_entries{:})

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
print(fig1, ['figs' filesep 'Figure1_deform_OM0' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_deform_OM0'], '-dsvg', '-r600');

%% Figure 2
fig2 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)+[width*100,0] width*100, height*100]); %<- Set size
    colormap(strainScheme)%     colorbar
        
axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15],'XTickLabel',{'','','',''} ,...
               'YLim',[-0.02,0.02]}; 
axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15],'XTickLabel',{'','','',''} ,...
               'YLim',[-1,1]*1e-3}; 
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15] ,...
               'YLim',[-1,1]*5e-5}; 
legend_entries = {'sphere','ellipsoid hor','ellipsoid ver'};

for j = [2,4,6]
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
        
    subplot(313); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end
    subplot(311)
        legend(legend_entries{:})
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
print(fig2, ['figs' filesep 'Figure1_deform_OM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure1_deform_OM10'], '-dsvg', '-r600');
