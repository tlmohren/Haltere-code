clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure3_strainData';
saveName = 'figure3_strainData';

renew_data_load = false
if renew_data_load
%     FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
%     FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(1).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om10';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om0';
    FEA(4).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om10';   
    FEA(5).name = 'Haltere_CraneFly_sphereVerOffset_Om0';
    FEA(6).name = 'Haltere_CraneFly_sphereVerOffset_Om10';   
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });        toc 
    end
    % Determine Circle locations
    for j = 1:length(FEA)
        circleDistance = 300;               % distance from base to haltere 
        circleRadius = 150;                 % radius of haltere   
        mindist =  min( abs( FEA(j).xyz(:,1) - circleDistance) );
        xMatch = find(  abs(FEA(j).xyz(:,1) - circleDistance) <= (mindist+1) );
        yMatch = find( round( abs( FEA(j).xyz(:,2) ), 7) == circleRadius );
        zMatch = find( round( abs( FEA(j).xyz(:,3) ), 7) == circleRadius );

        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end

%% deformation in angles 
fig1 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 2;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size

len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

strainLabel = {'$\epsilon$'};
axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2e-3,'YTick',[-1,0,1]*2e-3}; 
axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 
cols = linspecer(4); 

for j = [1,2]
    subplot(2,2,j); hold on 
        p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
        p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
        ylabel( strainLabel,'Rotation',0 );
        ax = gca();
        set(ax,axOptsStrainTop{:})
            set(p1,'Color',cols(1,:))
            set(p2,'Color',cols(2,:))
    subplot(2,2,j+2); hold on 
        p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
        p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
        xlabel('Time (s)')
        ylabel( strainLabel,'Rotation',0);
        ax = gca();
        set(ax,axOptsStrainSide{:})
            set(p3,'Color',cols(3,:))
            set(p4,'Color',cols(4,:))
end

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure3_strainOffsetHorOM0_OM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure3_strainOffsetHorOM0_OM10'], '-dsvg', '-r600');




%% deformation in angles 
fig2 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 2;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size

strainLabel = {'$\epsilon$'};
axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2e-3,'YTick',[-1,0,1]*2e-3}; 
axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 
cols = linspecer(4); 

for j = [3,4]
    subplot(2,2,j-2); hold on 
        p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
        p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
        ylabel( strainLabel,'Rotation',0 );
        ax = gca();
        set(ax,axOptsStrainTop{:})
            set(p1,'Color',cols(1,:))
            set(p2,'Color',cols(2,:))
    subplot(2,2,j); hold on 
        p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
        p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
        xlabel('Time (s)')
        ylabel( strainLabel,'Rotation',0);
        ax = gca();
        set(ax,axOptsStrainSide{:})
            set(p3,'Color',cols(3,:))
            set(p4,'Color',cols(4,:))
end

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure3_strainOffsetVerOM0_OM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure3_strainOffsetVerOM0_OM10'], '-dsvg', '-r600');



%% deformation in angles 
fig3 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 2;    % Height in inches
    set(fig3, 'Position', [fig3.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size

strainLabel = {'$\epsilon$'};
axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2e-3,'YTick',[-1,0,1]*2e-3}; 
axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 
cols = linspecer(4); 

for j = [5,6]
    subplot(2,2,j-4); hold on 
        p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
        p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
        ylabel( strainLabel,'Rotation',0 );
        ax = gca();
        set(ax,axOptsStrainTop{:})
            set(p1,'Color',cols(1,:))
            set(p2,'Color',cols(2,:))
    subplot(2,2,j-2); hold on 
        p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
        p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
        xlabel('Time (s)')
        ylabel( strainLabel,'Rotation',0);
        ax = gca();
        set(ax,axOptsStrainSide{:})
            set(p3,'Color',cols(3,:))
            set(p4,'Color',cols(4,:))
end

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig3,'InvertHardcopy','on');
set(fig3,'PaperUnits', 'inches');
papersize = get(fig3, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig3, 'PaperPosition', myfiguresize);
print(fig3, ['figs' filesep 'Figure3_strainOffsetVerSphereOM0_OM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig3, 'PaperPosition', myfiguresize);
print(fig3, ['figs' filesep 'Figure3_strainOffsetVerSphereOM0_OM10'], '-dsvg', '-r600');

% %%


% %%
% fig2 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size
% for j = [4]
%     subplot(211); hold on 
%         p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
%         p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
%         ylabel( strainLabel,'Rotation',0 );
%         ax = gca();
%         set(ax,axOptsStrainTop{:})
%     subplot(212); hold on 
%         p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
%         p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
%         xlabel('Time (s)')
%         ylabel( strainLabel,'Rotation',0);
%         ax = gca();
%         set(ax,axOptsStrainSide{:})
% end
% set(p1,'Color',cols(1,:))
% set(p2,'Color',cols(2,:))
% set(p3,'Color',cols(3,:))
% set(p4,'Color',cols(4,:))
% 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig2,'InvertHardcopy','on');
% set(fig2,'PaperUnits', 'inches');
% papersize = get(fig2, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig2, 'PaperPosition', myfiguresize);
% print(fig2, ['figs' filesep 'Figure3_strainOffsetHorOM10' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig2, 'PaperPosition', myfiguresize);
% print(fig2, ['figs' filesep 'Figure3_strainOffsetHorOM10'], '-dsvg', '-r600');

% %% 
% fig3 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig3, 'Position', [fig3.Position(1:2)-[width,height*1.5]*100 width*100, height*100]); %<- Set size
% 
% for j = [5,7]
%     subplot(211); hold on 
%         p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
%         p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
%         ylabel( strainLabel,'Rotation',0 );
%         ax = gca();
%         set(ax,axOptsStrainTop{:})
%     subplot(212); hold on 
%         p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
%         p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
%         xlabel('Time (s)')
%         ylabel( strainLabel,'Rotation',0);
%         ax = gca();
%         set(ax,axOptsStrainSide{:})
% end
% set(p1,'Color',cols(1,:))
% set(p2,'Color',cols(2,:))
% set(p3,'Color',cols(3,:))
% set(p4,'Color',cols(4,:))
% 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig3,'InvertHardcopy','on');
% set(fig3,'PaperUnits', 'inches');
% papersize = get(fig3, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig3, 'PaperPosition', myfiguresize);
% print(fig3, ['figs' filesep 'Figure3_strainOffsetVerOM0' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig3, 'PaperPosition', myfiguresize);
% print(fig3, ['figs' filesep 'Figure3_strainOffsetVerOM0'], '-dsvg', '-r600');
% 
% %%
% fig4 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig4, 'Position', [fig4.Position(1:2)-[0,height*1.5]*100 width*100, height*100]); %<- Set size
% 
% for j = [6,8]
%     subplot(211); hold on 
%         p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
%         p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
%         ylabel( strainLabel,'Rotation',0 );
%         ax = gca();
%         set(ax,axOptsStrainTop{:})
%     subplot(212); hold on 
%         p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
%         p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
%         xlabel('Time (s)')
%         ylabel( strainLabel,'Rotation',0);
%         ax = gca();
%         set(ax,axOptsStrainSide{:})
% end
% set(p1,'Color',cols(1,:))
% set(p2,'Color',cols(2,:))
% set(p3,'Color',cols(3,:))
% set(p4,'Color',cols(4,:))
% 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig4,'InvertHardcopy','on');
% set(fig4,'PaperUnits', 'inches');
% papersize = get(fig4, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig4, 'PaperPosition', myfiguresize);
% print(fig4, ['figs' filesep 'Figure3_strainOffsetVerOM10' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig4, 'PaperPosition', myfiguresize);
% print(fig4, ['figs' filesep 'Figure3_strainOffsetVerOM10'], '-dsvg', '-r600');
