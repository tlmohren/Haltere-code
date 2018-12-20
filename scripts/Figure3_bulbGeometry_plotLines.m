clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure3_bulbGeometry';
saveName = 'figure3_bulbGeometry';

% renew_data_load = true
renew_data_load = false
if renew_data_load
    FEA(1).name = 'Haltere_CraneFlyLowDensityWbulb_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_Sphere_Om10';
    FEA(3).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidVer_Om0';
    FEA(4).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_ellipsoidVer_Om10';
    FEA(5).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidHor_Om0';
    FEA(6).name = 'Haltere_CraneFlyLowDensityt8u7wBulb_ellipsoidHor_Om10';
    FEA(7).name =  'Haltere_CraneFlyLowDensitywBulb_SphereVerOffset_Om0';
    FEA(8).name =  'Haltere_CraneFlyLowDensitywBulb_SphereVerOffset_Om10';
    FEA(9).name =  'Haltere_CraneFlyLowDensitywBulb_ellipsoidVerOffset_Om0';
    FEA(10).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidVerOffset_Om10';   
    FEA(11).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidHorOffset_Om0';
    FEA(12).name = 'Haltere_CraneFlyLowDensitywBulb_ellipsoidHorOffset_Om10';

    for j =  1:length(FEA)
        tic
        [~, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
        [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc 
    end
    er = 0.1;
    circleDistance = 5000;               % distance from base to haltere 
    bulbA = [500,500, 353.543,353.543,    1000,1000,     500,500,   353.543,353.543,   1000,1000,   ] ; 
    bulbB = [500,500, 1000,1000,          353.543,353.54,500,500,   1000,1000,    353.543,353.543,] ; 
    yOffset = [0,0,0,0,0,0, 0,0, 0,0, 150, 150]; 
    zOffset = [0,0,0,0,0,0,  150, 150, 150, 150, 0,0]; 
    for j = 1:length(FEA)
        xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
        yMatch = find(   abs(abs( FEA(j).xyz(:,2) -yOffset(j)) - bulbA(j)   ) <er & ...
                     abs( abs( FEA(j).xyz(:,3)-zOffset(j) )   ) <er) ;
        zMatch = find(   abs(abs( FEA(j).xyz(:,3)-zOffset(j) ) - bulbB(j)   ) <er & ...
                     abs( abs( FEA(j).xyz(:,2)-yOffset(j) )   )<er );
        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    % Transpose deformation to haltere frame 
    for j =  1:length(FEA)
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
        FEA(j).xyzHaltereFrameMiddle = squeeze(  mean( FEA(j).xyzHaltereFrame( FEA(j).sideInds,:,:), 1)  )- ones(length(FEA(j).phi),1)*[0,yOffset(j),zOffset(j)];
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
fig3 = figure();
    width = 6;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig3, 'Position', [fig3.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
 

axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*1e-5};   
lineCols = linspecer(3); 
for k = 1:3
    for j = [1,2]
        jk = 2*(k-1)+j; 

        subplot(3,5, j ); hold on 
            p1 = plot( t_plot, FEA(jk).yAngle(It)  ,'color',lineCols(k,:) );
            ylabel( deformLabels{1} );
            ax = gca();
            set(ax,axOpts_dphi{:})

        subplot(3,5, j+5 ); hold on 
            p1 = plot( t_plot, FEA(jk).zAngle(It) ,'color',lineCols(k,:) );
            ylabel( deformLabels{2} );
            ax = gca();
            set(ax,axOpts_dtheta{:})

        subplot(3,5, j+ 10 ); hold on 
            p1 = plot( t_plot, FEA(jk).twistAngle(It)  ,'color',lineCols(k,:) );
            xlabel('Time (s)'); ylabel( deformLabels{3} );
            ax = gca();
            set(ax,axOpts_dgamma{:})
    end
end

           
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*5e-4};   

for k = 4:6
    
    for j = [1,2]
        jk = 2*(k-1)+j;
        subplot(3,5, j +2); hold on 
            p1 = plot( t_plot, FEA(jk).yAngle(It)  ,'color',lineCols(k-3,:) );
            ylabel( deformLabels{1} );
%             ax = gca();
%             set(ax,axOpts_dphi{:})

        subplot(3,5, j+7 ); hold on 
            p1 = plot( t_plot, FEA(jk).zAngle(It) ,'color',lineCols(k-3,:) );
            ylabel( deformLabels{2} );
            ax = gca();
            set(ax,axOpts_dtheta{:})

        subplot(3,5, j+12 ); hold on 
            p1 = plot( t_plot, FEA(jk).twistAngle(It)  ,'color',lineCols(k-3,:) );
            xlabel('Time (s)'); ylabel( deformLabels{3} );
            ax = gca();
            set(ax,axOpts_dgamma{:})
    end
end


axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*3e-5};   

for k = 4:6
    for j = [ 2]
        jk = 2*(k-1)+j;
        subplot(3,5, 15); hold on 
            p1 = plot( t_plot, FEA(jk).twistAngle(It)  ,'color',lineCols(k-3,:) );
            xlabel('Time (s)'); ylabel( deformLabels{3} );
            ax = gca();
            set(ax,axOpts_dgamma{:})
    end
end

%% Setting paper size for saving 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig3,'InvertHardcopy','on');
set(fig3,'PaperUnits', 'inches');
papersize = get(fig3, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig3, 'PaperPosition', myfiguresize);

print(fig3, ['figs' filesep 'Figure2_deformBulb' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig3, 'PaperPosition', myfiguresize);
print(fig3, ['figs' filesep 'Figure2_deformBulb'], '-dsvg', '-r600');


%% 
figure();
    subplot(311) 
    plot(  FEA(8).yAngle  -FEA(2).yAngle  )
    subplot(312)
    plot(  FEA(10).yAngle  -FEA(4).yAngle  )
    subplot(313)
    plot(  FEA(12).yAngle(1:150)  -FEA(6).yAngle(1:150)   )