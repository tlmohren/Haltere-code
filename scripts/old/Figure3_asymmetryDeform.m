clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure3_deform';
saveName = 'figure3_deform';

% renew_data_load = true
renew_data_load = false
if renew_data_load
%     FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
%     FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(1).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om10';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om0';
    FEA(4).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om10';   
    for j =  1:length(FEA)
        tic
        [~, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });
        [FEA(j).xyz, FEA(j).deform, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'u2','v2','w2'});
        [~, FEA(j).angles, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle','theta_angle'});
        toc 
    end
    % Determine Circle locations
    for j = 1:length(FEA)o
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
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size


deformLabels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};
legend_entries = {'ellipsoid hor', 'ellipsoid ver' };
len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;


axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-0.02,0.02]}; 
axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-1,1]*1e-3}; 
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*1e-3}; 
           
%            
% for j = [1,3]
%     subplot(3,1,1); hold on 
%         p1 = plot( t_plot, FEA(j).yAngle(It) );
%         ylabel( deformLabels{1} );
%         ax = gca();
%         set(ax,axOpts_dphi{:})
%         
%     subplot(312); hold on 
%         p1 = plot( t_plot, FEA(j).zAngle(It) );
%         ylabel( deformLabels{2} );
%         ax = gca();
%         set(ax,axOpts_dtheta{:})
%         
%     subplot(313); hold on 
%         p1 = plot( t_plot, FEA(j).twistAngle(It) );
%         xlabel('Time (s)'); ylabel( deformLabels{3} );
%         ax = gca();
%         set(ax,axOpts_dgamma{:})
% end
%

% figure()
for j = [1,2]
    subplot(3,2, j ); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(3,2, j+2 ); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(3,2, j+4 ); hold on 
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
print(fig1, ['figs' filesep 'Figure3_deformOffsetHorOM0_OM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure3_deformOffsetHorOM0_OM10'], '-dsvg', '-r600');

%% deformation in angles 
fig2 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size

% for j = [2,6]
%     subplot(3,1,1); hold on 
%         p1 = plot( t_plot, FEA(j).yAngle(It) );
%         ylabel( deformLabels{1} );
%         ax = gca();
%         set(ax,axOpts_dphi{:})
%         
%     subplot(312); hold on 
%         p1 = plot( t_plot, FEA(j).zAngle(It) );
%         ylabel( deformLabels{2} );
%         ax = gca();
%         set(ax,axOpts_dtheta{:})
%         
%     subplot(313); hold on 
%         p1 = plot( t_plot, FEA(j).twistAngle(It) );
%         xlabel('Time (s)'); ylabel( deformLabels{3} );
%         ax = gca();
%         set(ax,axOpts_dgamma{:})
% end

axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-0.02,0.02]}; 
axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
               'YLim',[-1,1]*1e-3}; 
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*1e-4}; 

for j = 3:4
    subplot(3,2, j-2 ); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(3,2, j ); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(3,2, j+2 ); hold on 
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

print(fig2, ['figs' filesep 'Figure3_deformOffsetVerOM0_OM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure3_deformOffsetVerOM0_OM10'], '-dsvg', '-r600');


%% deformation in angles 
% fig3 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 3;    % Height in inches
%     set(fig3, 'Position', [fig3.Position(1:2)-[width,height*1.5]*100 width*100, height*100]); %<- Set size
% 
% 
% deformLabels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};
% legend_entries = {'ellipsoid hor', 'ellipsoid ver' };
% len = 101;
% start = 35;
% It = start:(start+len-1);
% t_plot = (0:len-1)*0.001;
% 
% 
% axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
%                'YLim',[-0.02,0.02]}; 
% axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,...
%                'YLim',[-1,1]*1e-3}; 
% axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
%                'YLim',[-1,1]*1e-3}; 
%            
%            
% for j = [1,5,3]
%     subplot(3,1,1); hold on 
%         p1 = plot( t_plot, FEA(j).yAngle(It) );
%         ylabel( deformLabels{1} );
%         ax = gca();
%         set(ax,axOpts_dphi{:})
%         
%     subplot(312); hold on 
%         p1 = plot( t_plot, FEA(j).zAngle(It) );
%         ylabel( deformLabels{2} );
%         ax = gca();
%         set(ax,axOpts_dtheta{:})
%         
%     subplot(313); hold on 
%         p1 = plot( t_plot, FEA(j).twistAngle(It) );
%         xlabel('Time (s)'); ylabel( deformLabels{3} );
%         ax = gca();
%         set(ax,axOpts_dgamma{:})
% end
% 
% %% Setting paper size for saving 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig3,'InvertHardcopy','on');
% set(fig3,'PaperUnits', 'inches');
% papersize = get(fig3, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig3, 'PaperPosition', myfiguresize);
% print(fig3, ['figs' filesep 'Figure3_deformOffsetHorOM0' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig3, 'PaperPosition', myfiguresize);
% print(fig3, ['figs' filesep 'Figure3_deformOffsetHorOM0'], '-dsvg', '-r600');
% 
% 
% %% deformation in angles 
% fig4 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 3;    % Height in inches
%     set(fig4, 'Position', [fig4.Position(1:2)-[0,height*1.5]*100 width*100, height*100]); %<- Set size
% 
% for j = [2,6,4]
%     subplot(3,1,1); hold on 
%         p1 = plot( t_plot, FEA(j).yAngle(It) );
%         ylabel( deformLabels{1} );
%         ax = gca();
%         set(ax,axOpts_dphi{:})
%         
%     subplot(312); hold on 
%         p1 = plot( t_plot, FEA(j).zAngle(It) );
%         ylabel( deformLabels{2} );
%         ax = gca();
%         set(ax,axOpts_dtheta{:})
%         
%     subplot(313); hold on 
%         p1 = plot( t_plot, FEA(j).twistAngle(It) );
%         xlabel('Time (s)'); ylabel( deformLabels{3} );
%         ax = gca();
%         set(ax,axOpts_dgamma{:})
% end
% 
% %% Setting paper size for saving 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig4,'InvertHardcopy','on');
% set(fig4,'PaperUnits', 'inches');
% papersize = get(fig4, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig4, 'PaperPosition', myfiguresize);
% print(fig4, ['figs' filesep 'Figure3_deformOffsetHorOM10' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig4, 'PaperPosition', myfiguresize);
% print(fig4, ['figs' filesep 'Figure3_deformOffsetHorOM10'], '-dsvg', '-r600');
% 
% 
% 
% %% appendix figure 
% j = 1;
% 
% fig2 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig2, 'Position', [fig2.Position(1:2)+[width*100,height*100] width*100, height*100]); %<- Set size
%         plot(t_plot,FEA(j*2-1).twistAngle(It) - FEA(j*2).twistAngle(It) )
% %         plot(t(It),FEA(j*2).twistAngle(It) , lineSpec{j*2})
%         xlabel('Time (s)'); ylabel( deformLabels{3} );
%         ax = gca();
%         set(ax,axOpts{:})