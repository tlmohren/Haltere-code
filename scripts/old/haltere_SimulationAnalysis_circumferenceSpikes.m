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

% Ind = circleIndices(I_sort);
circlePoints = circleIndices(I_sort); 
sidePoints = circlePoints( find( mod(V,90) == 0));

%% 
% pre plot decisions 
width = 6;     % Width in inches,   find column width in paper 
height = 5;    % Height in inches
fsz = 8;      % Fontsize

fsz = 10;
set(0,'DefaultAxesFontSize',fsz)% .
set(0,'DefaultLegendFontSize',fsz)% .

fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

[tI,nS,s ] = size( sim.Om0.strain);
t = ( (1:tI) - tI+50) /25;
col = linspecer(2);

%% 
%% Create neural encoding functions      
STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );

%% 

[tL, n, nS] = size( sim.Om0.strain(:, sidePoints,1) );   
fSamp = 1000;
subSamp =10;
STAt = linspace(-39,0,40*subSamp);

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
% calib_param = [0.0030    0.0542    0.0029    0.0515];
calib_param = [0.0030    0.0204    0.0377    0.0497    0.0542    0.0505    0.0393    0.0223    0.0029    0.0211    0.0372    0.0479    0.0515    0.0472    0.0360    0.0197];
%% apply neural encoding 

for j = 1:length(circlePoints)
    strainTemp = sim.Om0.strain(:,circlePoints(j),1);
    t = (1:tL)/fSamp;
    tNew = linspace(t(1),t(end), tL*subSamp ) ; 
    strainInterp(:,j) =  interp1(t,strainTemp ,tNew,'spline');

    STA = STAfun(STAt);
    strainConv(:,j) = conv( [zeros(1,length(STA)-1),strainInterp(:,j)'], fliplr( STA), 'valid');
%     calib_param(j) = max(  strainConv(:,j) );
    sim.Om0.pFire(:,j) = NLDfun( strainConv(:,j)/calib_param(j) );
    sim.Om0.spikeInds{j} = findSpikes( sim.Om0.pFire(:,j) ); 
end


for j = 1:length(circlePoints)
    strainTemp = sim.Om10.strain(:,circlePoints(j),1);
    t = (1:tL)/fSamp;
    tNew = linspace(t(1),t(end), tL*subSamp ) ; 
    strainInterp(:,j) =  interp1(t,strainTemp ,tNew,'spline');

    STA = STAfun(STAt);
    strainConv(:,j) = conv( [zeros(1,length(STA)-1),strainInterp(:,j)'], fliplr( STA), 'valid');
%     calib_param(j) = max(  strainConv(:,j) );
    sim.Om10.pFire(:,j) = NLDfun( strainConv(:,j)/calib_param(j) );
    sim.Om10.spikeInds{j} = findSpikes( sim.Om10.pFire(:,j) ); 
    
    
end
%% 
for j = 1:length(circlePoints)
   dI = sim.Om0.spikeInds{j}(end) - sim.Om10.spikeInds{j}(end);
   dT(j) = dI/1e4;
end

%% 

axCircle= {
    'XLim',[-0.5,2],...
%     'XTick',[0:1:2],...
}; 
    

subNr = [1:2:15,2:2:16];
subNr = reshape(1:27,3,9)';


tS = ( (1:tI) - tI+50) /25;
tNon = linspace(t(1),t(end), tL*subSamp )*40 -1.24 ; 

fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size

for j = 6:12
    subplot(9,3, subNr(j-4,1)); hold on
    
    plot( tS, sim.Om0.strain(:, circlePoints(j),1) ,'Color',col(1,:) )
    plot( tS, sim.Om10.strain(:, circlePoints(j),1) ,'Color',col(2,:) )
    axis off
    lim = max( sim.Om10.strain(:, circlePoints(j),1)); 
    locs = tNon(sim.Om10.spikeInds{j});
    plot( [1,1]'*locs   ,  [lim*1, lim*0.3]'*[1,1] ,'r','LineWidth',1)
    try
%         plot( tNon(sim.Om0.spikeInds{j}) , zeros(length(sim.Om0.spikeInds{j}),1)','or')
        locs = tNon(sim.Om0.spikeInds{j});
        plot( [1,1]'*locs   ,  [-lim*0.3, -lim*1]'*[1,1] ,'b','LineWidth',1)
    end
    plot( [-0.5,2],[0,0],'Color',[1,1,1]*0.5)
    text(2,0,num2str( dT(j)*1e3 ) )
    
    ax = gca();
    set(ax,axCircle{:})
end

for j = 5
    subplot(9,3, subNr(j-4,2)); hold on
    
    plot( tS, sim.Om0.strain(:, circlePoints(j),1) ,'Color',col(1,:) )
    plot( tS, sim.Om10.strain(:, circlePoints(j),1) ,'Color',col(2,:) )
    axis off
    lim = max( sim.Om10.strain(:, circlePoints(j),1)); 
    locs = tNon(sim.Om10.spikeInds{j});
    plot( [1,1]'*locs   ,  [lim*1, lim*0.3]'*[1,1] ,'r','LineWidth',1)
    try
%         plot( tNon(sim.Om0.spikeInds{j}) , zeros(length(sim.Om0.spikeInds{j}),1)','or')
        locs = tNon(sim.Om0.spikeInds{j}(end));
        plot( [1,1]'*locs   ,  [-lim*0.3, -lim*1]'*[1,1] ,'b','LineWidth',1)
    end
    plot( [-0.5,2],[0,0],'Color',[1,1,1]*0.5)
    text(2,0,num2str( dT(j)*1e3 ) )
    ax = gca();
    set(ax,axCircle{:})
end


for j = 13
    subplot(9,3, subNr(j-4,2)); hold on

    
    plot( tS, sim.Om0.strain(:, circlePoints(j),1) ,'Color',col(1,:) )
    plot( tS, sim.Om10.strain(:, circlePoints(j),1) ,'Color',col(2,:) )
    axis off
    lim = max( sim.Om10.strain(:, circlePoints(j),1)); 
    locs = tNon(sim.Om10.spikeInds{j});
    plot( [1,1]'*locs   ,  [lim*1, lim*0.3]'*[1,1] ,'r','LineWidth',1)
    try
%         plot( tNon(sim.Om0.spikeInds{j}) , zeros(length(sim.Om0.spikeInds{j}),1)','or')
        locs = tNon(sim.Om0.spikeInds{j}(end));
        plot( [1,1]'*locs   ,  [-lim*0.3, -lim*1]'*[1,1] ,'b','LineWidth',1)
    end
    plot( [-0.5,2],[0,0],'Color',[1,1,1]*0.5)
    text(2,0,num2str( dT(j)*1e3 ) )
    ax = gca();
    set(ax,axCircle{:})
end



count = 2;
for j = fliplr([14:16,1:4])
    subplot(9,3, subNr(count,3)); hold on
    count = count+1;
    
    plot( tS, sim.Om0.strain(:, circlePoints(j),1) ,'Color',col(1,:) )
    plot( tS, sim.Om10.strain(:, circlePoints(j),1) ,'Color',col(2,:) )
    
    lim = max( sim.Om10.strain(:, circlePoints(j),1)); 
    locs = tNon(sim.Om10.spikeInds{j});
    plot( [1,1]'*locs   ,  [lim*1, lim*0.3]'*[1,1] ,'r','LineWidth',1)
    try
%         plot( tNon(sim.Om0.spikeInds{j}) , zeros(length(sim.Om0.spikeInds{j}),1)','or')
        locs = tNon(sim.Om0.spikeInds{j}(end));
        plot( [1,1]'*locs   ,  [-lim*0.3, -lim*1]'*[1,1] ,'b','LineWidth',1)
    end
    plot( [-0.5,2],[0,0],'Color',[1,1,1]*0.5)
    
    text(2,0,num2str( dT(j)*1e3 ) )
    axis off
    grid on
    ax = gca();
    set(ax,axCircle{:})
end



axCirc = {...
    'PlotBoxAspectRatio',[3 4 4],...
    'DataAspectRatio',[1 1 1],...
    'XLim',[-1,1],...
    'YLim',[-1,1],...
    };
subplot(9,3,subNr(4:6,2) ); hold on
    x = sin(  linspace(0,1,100)*2*pi );
    y = cos( linspace(0,1,100)*2*pi );
    plot(x,y,'k')

    x = sin(  linspace(0,1,17)*2*pi );
    y = cos( linspace(0,1,17)*2*pi );
    scatter(x,y,'ro','filled')
    ax = gca();
    set(ax,axCirc{:} )
    axis off

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
print(fig1, ['figs' filesep 'Figure_circleSpikes' ], '-dpng', '-r600');

% total hack, why does saving to svg scale image up???
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);

print(fig1, ['figs' filesep 'Figure_circleSpikes' ], '-dsvg');
