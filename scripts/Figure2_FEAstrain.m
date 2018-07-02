clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure2_strainData';
saveName = 'figure2_strainData';

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

strainLabel = {'$\epsilon$'};
axOptsStrainTop = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','','',''},'YLim',[-1,1]*2e-3,'YTick',[-1,0,1]*2e-3}; 
axOptsStrainSide = {'XGrid','On','XLim',[0,0.15],'XTick',[0:0.05:0.15],'YLim',[-1,1]*1e-4}; 

cols = linspecer(4); 
for j = [1,3,5]
    subplot(211); hold on 
        p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It) );
        p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
        ylabel( strainLabel,'Rotation',0 );
        ax = gca();
        set(ax,axOptsStrainTop{:})
    subplot(212); hold on 
        p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
        p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
        xlabel('Time (s)')
        ylabel( strainLabel,'Rotation',0);
        ax = gca();
        set(ax,axOptsStrainSide{:})
end
set(p1,'Color',cols(1,:))
set(p2,'Color',cols(2,:))
set(p3,'Color',cols(3,:))
set(p4,'Color',cols(4,:))
        
%% Setting paper size for saving 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_strainOM0' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_strainOM0'], '-dsvg', '-r600');

%% deformation in angles 
fig2 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 2;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size

for j = [2,4,6]
    subplot(211); hold on 
        p1 = plot(t_plot, FEA(j).strain( FEA(j).topInds(1), It ) );
        p2 = plot(t_plot, FEA(j).strain( FEA(j).topInds(2), It) );
        ylabel( strainLabel ,'Rotation',0);
        ax = gca();
        set(ax,axOptsStrainTop{:})
    subplot(212); hold on 
        p3 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(1), It) );
        p4 = plot(t_plot, FEA(j).strain( FEA(j).sideInds(2), It) );
        ylabel( strainLabel,'Rotation',0 );
        xlabel('Time (s)')
        ax = gca();
        set(ax,axOptsStrainSide{:})
end
set(p1,'Color',cols(1,:))
set(p2,'Color',cols(2,:))
set(p3,'Color',cols(3,:))
set(p4,'Color',cols(4,:))       
        
%% Setting paper size for saving 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure2_strainOM10' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure2_strainOM10'], '-dsvg', '-r600');