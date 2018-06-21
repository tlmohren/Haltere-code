%-------------------------------
% TMohren
% load comsol data and plot resulting strains
% 2017-08-09
%------------------------------
clc; clear all; close all
addpathFolderStructureHaltere()

sim.name = 'Haltere_CraneFly_Sphere';
load(['data' filesep sim.name filesep sim.name '_allData'])

%% Simulation parameters
circleDistance = 300;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere      
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

%% 
%  % plot strain at pointIndices
%% Create neural encoding functions      
STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );


    
%% 

% [tL, n, nS] = size( sim.Om0.strain  )
[tL, n, nS] = size( sim.Om0.strain(:, sidePoints,1) )
    
fSamp = 1000;
subSamp =10;
STAt = linspace(-39,0,40*subSamp);

parse_extra = 5;
% for j = 1:length(simName)

%% 


axConv= {
    'Xgrid','on',...
    'Ygrid','on',...
    'box','off',...
    'XLim',[-0.5,2],...
    'XTick',[0:1:2],...
%     'XTicklabel',{'-40','','-20','','0'},...
%     'YLim',[-1.2,0.7],...
%     'YTick',[],...
    };

calib_param = [0.0030    0.0542    0.0029    0.0515];
%% apply neural encoding 


for j = 1:4
    strainTemp = sim.Om0.strain(:,sidePoints(j),1);
    t = (1:tL)/fSamp;
    tNew = linspace(t(1),t(end), tL*subSamp ) ; 
    strainInterp(:,j) =  interp1(t,strainTemp ,tNew,'spline');


    STA = STAfun(STAt);
    strainConv(:,j) = conv( [zeros(1,length(STA)-1),strainInterp(:,j)'], fliplr( STA), 'valid');
    calib_param(j) = max(  strainConv(:,j) );
    pFire(:,j) = NLDfun( strainConv(:,j)/calib_param(j) );
    
    spikeInds{j} = findSpikes( pFire(:,j) ); 
end

% for j = 1:4
%     spikeInds{j} = findSpikes( pFire(:,j) ); 
% end
% find spikes 
%%
j =2 ;
tNon = linspace(t(1),t(end), tL*subSamp )*40 -1.24 ; 

figure();
subplot(311)
    plot(tNon,strainInterp(:,j))
    ax = gca();
    set(ax,axConv{:})
subplot(312)
    plot(tNon,strainConv(:,j))
    ax = gca();
    set(ax,axConv{:})
subplot(313) 
    plot(tNon,pFire(:,j));hold on
    plot( (tNon(spikeInds{j})'*[1,1])' , [0,1]'*ones(length(spikeInds{j}),1)','k','LineWidth',1 )
    plot( tNon(spikeInds{j}) , ones(length(spikeInds{j}),1)','ok')
    ax = gca();
    set(ax,axConv{:})
    
%% 
tS = ( (1:tI) - tI+50) /25;

fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
subplot(411); hold on
    plot( tS, sim.Om0.strain(:, sidePoints(2),1) ,'Color',col(1,:) )
    ax = gca();
    set(ax,axOpts{:})
    ylabel('Top Strain','Rotation',0); 
    sim.Om0.strain(:, sidePoints(2),1) 
     plot( tNon(spikeInds{2}) , zeros(length(spikeInds{2}),1)','ok')
subplot(412); hold on
    plot( tS, sim.Om0.strain(:, sidePoints(4),1),'Color',col(1,:) )
    ax = gca();
    set(ax,axOpts{:})
    ylabel('Bottom Strain','Rotation',0)
     plot( tNon(spikeInds{4}) , zeros(length(spikeInds{4}),1)','ok')
     
     
subplot(413); hold on
    plot( tS, sim.Om0.strain(:, sidePoints(1),1),'Color',col(2,:) )
    ax = gca();
    set(ax,axOpts{:})
    ylabel('Left Strain','Rotation',0)
     plot( tNon(spikeInds{1}) , zeros(length(spikeInds{1}),1)','ok')
    
subplot(414); hold on
    plot( tS, sim.Om0.strain(:, sidePoints(3),1),'Color',col(2,:) )
    ax = gca();
    set(ax,axOpts{:})
    xlabel('Haltere Beats (t/T)'); ylabel('Right Strain','Rotation',0)
     plot( tNon(spikeInds{3}) , zeros(length(spikeInds{3}),1)','ok')
%     
%     
%     
%     
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
print(fig1, ['figs' filesep 'Figure_spikesF' ], '-dpng', '-r600');

% total hack, why does saving to svg scale image up???
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);

print(fig1, ['figs' filesep 'Figure_spikesF' ], '-dsvg');
