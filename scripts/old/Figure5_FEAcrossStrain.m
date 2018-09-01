clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure5_strainData';
saveName = 'figure5_strainData';

renew_data_load = false 
% renew_data_load = true
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om0';
    FEA(2).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om10';
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX' });        toc 
    end
    sideL = 243.435; 
    circleDistance = 300;               % distance from base to haltere 
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
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end
% 
%% deformation in angles 
% 
% fig1 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
% len = 101;
% start = 35;
% It = start:(start+len-1);
% t_plot = (0:len-1)*0.001;
% 
% axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2.5e-3,'YTick',[-1,0,1]*2e-3}; 
% axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 
% 
% for j = [1,2]
%     subplot(2,2,j); hold on 
%         plot(t_plot, FEA(j).strain( FEA(j).topInds, It) )
%         ax = gca();
%         set(ax,axOptsStrainTop{:})
%         ylabel('$\epsilon$','rotation',0)
%     subplot(2,2,j+2); hold on 
%         plot(t_plot, FEA(j).strain( FEA(j).sideInds, It) )
%         ax = gca();
%         set(ax,axOptsStrainSide{:})
%         xlabel('Time (s)')
%         ylabel('$\epsilon$','rotation',0)
% end
% 
% %% Setting paper size for saving 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig1,'InvertHardcopy','on');
% set(fig1,'PaperUnits', 'inches');
% papersize = get(fig1, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig1, 'PaperPosition', myfiguresize);
% print(fig1, ['figs' filesep 'Figure5_strainPlotOM0' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig1, 'PaperPosition', myfiguresize);
% print(fig1, ['figs' filesep 'Figure5_strainPlotOM0'], '-dsvg', '-r600');
% % 




%% deformation in angles 

fig2 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 2;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size
len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

cols = linspecer(4); 
axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2.5e-3,'YTick',[-1,0,1]*2e-3}; 
axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 

for j = [1,2]
    subplot(2,2,j); hold on 
        p1= plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
        p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
        ax = gca();
        set(ax,axOptsStrainTop{:})
        ylabel('$\epsilon$','rotation',0)
            set(p1,'Color',cols(1,:))
            set(p2,'Color',cols(2,:))
    subplot(2,2,j+2); hold on 
        p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
        p4 =plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
        ax = gca();
        set(ax,axOptsStrainSide{:})
        xlabel('Time (s)')
        ylabel('$\epsilon$','rotation',0)
            set(p3,'Color',cols(3,:))
            set(p4,'Color',cols(4,:))
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
print(fig2, ['figs' filesep 'Figure5_strainPlotOM0_OM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure5_strainPlotOM0_OM10'], '-dsvg', '-r600');
% 