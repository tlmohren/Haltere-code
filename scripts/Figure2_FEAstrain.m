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


% 
% 
% 
% 
%% deformation in angles 

fig1 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

t = 0:0.001:0.35;
labels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};
axOpts1 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 
axOpts2 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 
axOpts3 = {'XGrid','On','XLim',[0,0.2],'XTick',[0:0.05:0.2]}; 

lineSpec = {'-','o','--','+','-.','x'};

legend_entries = {'sphere','sphere 10','ellipsoid hor','ellipsoid hor 10', 'ellipsoid ver', 'ellipsoid ver 10'};

for j = 1:length(FEA)/2
    subplot(211); hold on 
        plot(t, FEA(j*2-1).strain( FEA(j*2-1).topInds, :) , lineSpec{j*2-1})
        plot(t, FEA(j*2).strain( FEA(j*2-1).topInds, :) , lineSpec{j*2-1})
    j
    subplot(212); hold on 
        plot(t, FEA(j*2-1).strain( FEA(j*2-1).sideInds, :) , lineSpec{j*2-1})
        plot(t, FEA(j*2).strain( FEA(j*2-1).sideInds, :) , lineSpec{j*2-1})
end
% 
%% Setting paper size for saving 

set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_strainPlot' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_strainPlot'], '-dsvg', '-r600');
% 