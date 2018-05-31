%-------------------------------
% TMohren
% load comsol data and plot resulting strains
% 2017-08-09
%------------------------------
clc; clear all; close all
addpathFolderStructureHaltere()

sim.name = 'Haltere_CraneFly_Sphere';
load(['data' filesep sim.name filesep sim.name '_allData'])

% Simulation parameters
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere
strainPoints = [circleDistance, 0,   circleRadius; ...
                circleDistance, 0,   - circleRadius; ...
               circleDistance,  circleRadius, 0;...
               circleDistance, -circleRadius, 0];




%% Analyze strain
circleIndices = [];

% pointIndices= findPointIndices( round(sim.strainXYZ,7) , strainPoints );
circleIndices= findCircleIndices( round(sim.strainXYZ,7) , circleDistance,circleRadius);
xyz = sim.strainXYZ;
angle = atan2( xyz(3,circleIndices), xyz(2,circleIndices));
angleDeg = rad2deg(angle)-180;
angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
[V,I_sort] = sort(angleDeg,'ascend');

Ind = circleIndices(I_sort);
sidePoints = Ind( find( mod(V,90) == 0));


%% 

dotStyle = {'Marker','.'};
markerStyle = {'Marker','o','MarkerFaceColor','red','CData',[1,0,0]};
    

% pre plot decisions 
width = 3.5;     % Width in inches,   find column width in paper 
height = 4;    % Height in inches
fsz = 8;      % Fontsize

fsz = 10;
           
set(0,'DefaultAxesFontSize',fsz)% .
set(0,'DefaultLegendFontSize',fsz)% .

fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

[tI,nS,s ] = size( sim.Om0.strain);
t = ( (1:tI) - tI+50) /25;
% axisOptsFig1A = {'xtick',[0:10:30,errLocFig1A  ],'xticklabel',{'0','10','20','30', '\bf 1326'},...
%     'ytick',0.5:0.25:1 ,'xlim', [0,errLocFig1A+2],'ylim',[0.4,1] ,...
%     'LabelFontSizeMultiplier',1};
axOpts = {'XLim',[-0.5,2],...
    'XTick',[0,1,2],...
    'Xgrid','on',...
    'Ygrid','on',...
    };

col = linspecer(2);
% col = {[1,1,1]*100/255,'-r'};
% dotcol = {'.k','.r'}; 
% % hold on
%% 
%  % plot strain at pointIndices
subplot(411); hold on
    plot( t, sim.Om0.strain(:, sidePoints(2),1) ,'Color',col(1,:) )
    ax = gca();
    set(ax,axOpts{:})
    ylabel('Top Strain','Rotation',0); 
subplot(412); hold on
    plot( t, sim.Om0.strain(:, sidePoints(4),1),'Color',col(1,:) )
    ax = gca();
    set(ax,axOpts{:})
    ylabel('Bottom Strain','Rotation',0)
    
    
subplot(413); hold on
    plot( t, sim.Om0.strain(:, sidePoints(1),1),'Color',col(2,:) )
    ax = gca();
    set(ax,axOpts{:})
    ylabel('Left Strain','Rotation',0)
    
subplot(414); hold on
    plot( t, sim.Om0.strain(:, sidePoints(1),1),'Color',col(2,:) )
    ax = gca();
    set(ax,axOpts{:})
    xlabel('Haltere Beats (t/T)'); ylabel('Right Strain','Rotation',0)


%% Setting paper size for saving 

set(0,'DefaultAxesFontSize',8)% .
set(0,'DefaultLegendFontSize',8)% .


set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
tightfig;
% % Here we preserve the size of the image when we save it.


set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);

% Saving figure 
print(fig1, ['figs' filesep 'Figure_Strains' ], '-dpng', '-r600');

% total hack, why does saving to svg scale image up???
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);

print(fig1, ['figs' filesep 'Figure_Strains' ], '-dsvg');
