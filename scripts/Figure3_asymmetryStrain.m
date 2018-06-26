clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure3_strainData';
saveName = 'figure3_strainData';

renew_data_load = false
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om0';
    FEA(4).name = 'Haltere_CraneFly_ellipsoidHorOffset_Om10';
    FEA(5).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om0';
    FEA(6).name = 'Haltere_CraneFly_ellipsoidVerOffset_Om10';   
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


% 
% 
% 
% 
%% deformation in angles 

fig1 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

t = 0:0.001:0.35;


labels = {'Top and bottom strain','left and right strain'};
axOpts1 = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.2]}; 
axOpts2 = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.2]}; 
axOpts3 = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.2]}; 

lineSpec = {'o','-','+','-','x','-'};

legend_entries = {'sphere0 top','sphere0 bottom',...
                'sphere10 top','sphere10 bottom',...
                'ellipsoidHorOffset0 top','ellipsoidHorOffset0 bottom',...
                'ellipsoidHorOffset10 top','ellipsoidHorOffset10 bottom',...
                'ellipsoidVerOffset0 top','ellipsoidVerOffset0 bottom',...
                'ellipsoidVerOffset10 top','ellipsoidVerOffset10 bottom',...
                };
It = 1:160;


for j = 1:length(FEA)/2
    subplot(211); hold on 
        plot(t(It), FEA(j*2-1).strain( FEA(j*2).topInds, (It)) , lineSpec{j*2-1})
        plot(t(It), FEA(j*2).strain( FEA(j*2).topInds, (It)) , lineSpec{j*2})
        ylabel( labels{2} );
        ax = gca();
        set(ax,axOpts2{:})
    subplot(212); hold on 
        plot(t(It), FEA(j*2-1).strain( FEA(j*2).sideInds, (It)) , lineSpec{j*2-1})
        plot(t(It), FEA(j*2).strain( FEA(j*2).sideInds, (It)) , lineSpec{j*2})
        ylabel( labels{2} );
        ax = gca();
        set(ax,axOpts2{:})
end
legend(legend_entries)
        
        
%% deformation in angles 

fig2 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size

t = 0:0.001:0.35;


labels = {'Top and bottom strain','left and right strain'};
axOpts1 = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.2]}; 
axOpts2 = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.2]}; 
axOpts3 = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.2]}; 

lineSpec = {'o','-','+','-','x','-'};

legend_entries = {'sphere0 top',...
                'sphere10 top',...
                'ellipsoidHorOffset0 bottom',...
                'ellipsoidHorOffset10 bottom',...
               'ellipsoidVerOffset0 bottom',...
                'ellipsoidVerOffset10 bottom',...
                };
It = 1:160;


for j = 1:length(FEA)/2
    subplot(211); hold on 
        plot(t(It), -FEA(j*2-1).strain( FEA(j*2).topInds(1), (It) ) +FEA(j*2-1).strain( FEA(j*2).topInds(2) , (It)) , lineSpec{j*2-1})
        plot(t(It), -FEA(j*2).strain( FEA(j*2).topInds(1), (It))+ FEA(j*2).strain( FEA(j*2).topInds(2), (It)) , lineSpec{j*2})
        ylabel( labels{2} );
        ax = gca();
        set(ax,axOpts2{:})
    subplot(212); hold on 
        plot(t(It), FEA(j*2-1).strain( FEA(j*2).sideInds(1), (It)) -  FEA(j*2-1).strain( FEA(j*2).sideInds(2), (It)) , lineSpec{j*2-1})
        plot(t(It), FEA(j*2).strain( FEA(j*2).sideInds(1), (It))- FEA(j*2).strain( FEA(j*2).sideInds(2), (It))  , lineSpec{j*2})
        ylabel( labels{2} );
        ax = gca();
        set(ax,axOpts2{:})
end
legend(legend_entries)
%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure3_strainOffsetPlot' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure3_strainOffsetPlot'], '-dsvg', '-r600');
% 