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

% pre plot decisions 
width = 3.5;     % Width in inches,   find column width in paper 
height = 4;    % Height in inches
fsz = 10;
           
set(0,'DefaultAxesFontSize',fsz)% .
set(0,'DefaultLegendFontSize',fsz)% .



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
% 

axSTA= {
    'Xgrid','on',...
    'Ygrid','on',...
    'XTick',[-40:10:0],...
    'XTicklabel',{'-40','','-20','','0'},...
    'YLim',[-1.2,0.7],...
    'YTick',[],...
    };

col = linspecer(2);


%% Create neural encoding functions      
STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );


%% 

fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

% 

axSTA= {
      'Visible','off',...
%       'box', 'off',...
%       'XTick',[],...
%       'YTick',[],...
%     'Xgrid','on',...
%     'Ygrid','on',...
%     'XTick',[-40:20:0],...
%     'YLim',[-1.2,0.7],...
%     'YTick',-1:0.5:1,...
%     'XTicklabel',{'-40','','-20','','0'},...
    };
axNLD= {
      'Visible','off',...
%     'Xgrid','on',...
%     'Ygrid','on',...
%     'XTick',[-1:1:1],...
%     'YLim',[-0.1,1.1],...
%     'YTick',0:0.5:1,...
%     'XTicklabel',{'-40','','-20','','0'},...
    };

subplot(211)
    t = -40:0.1:0;
    plot(t,STAfun(t),'k','LineWidth', 2) 
%     ylabel('Strain','Rotation',0)
    ax = gca();
    set(ax,axSTA{:})
%     axis off 
    
subplot(212)
    s = -1:0.01:1;
    plot(s,NLDfun(s),'k','LineWidth', 2)
%     ylabel('Prob. Firing','Rotation',0)
    ax = gca();
    set(ax,axNLD{:})
%     axis off 
    
    
%% 

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
print(fig1, ['figs' filesep 'Figure_Encoder' ], '-dpng', '-r600');

% total hack, why does saving to svg scale image up???
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);

print(fig1, ['figs' filesep 'Figure_Encoder' ], '-dsvg');