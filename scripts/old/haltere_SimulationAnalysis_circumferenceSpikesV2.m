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
   dT(j) = abs(dI/1e4 );
end

%% 




pY = 9;
pX = 3; 
% subNr = [1:2:15,2:2:16];
subNr = reshape(1: pY*pX,pX, pY)';


tS = ( (1:tL) - tL+50) /25;
tNon = linspace(t(1),t(end), tL*subSamp )*40 -1.24 ; 


col = linspecer(2);


fig1 = figure();
set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size




axCircle= {
    'XLim',[1,2],...
    'XTick',0:1:2,...
    'Visible','off',...
    'YLim',[-1,1],...
    'Ygrid','on',...
}; 
    
lw= 2; 
slw = 2;
blw = 1;
%
for j = 6:12
    subplot( pY, pX, subNr(j-4,1)); hold on
    
    plot( tNon, sin(2*pi*tNon),':','Color', [1,1,1]*0.1)
    
    plot( tNon, -sim.Om0.pFire(:, j,1) ,'k' ,'LineWidth',lw)
    plot(  tNon, sim.Om10.pFire(:, j,1) ,'Color',col(2,:) ,'LineWidth',lw)
%     axis off
    lim = max( sim.Om10.pFire(:, j ,1)); 
    locs = tNon(sim.Om10.spikeInds{j});
    plot( [1,1]'*locs   ,  [lim, 0]'*0.6,'k','LineWidth',slw)
    try
        locs = tNon(sim.Om0.spikeInds{j});
        plot( [1,1]'*locs   ,  [0, -lim]'*0.6 ,'k','LineWidth',slw)
    end
    plot( [-0.5,2],[0,0],'Color',[1,1,1]*0.5 )
%     text(2.1,0,num2str( dT(j)*1e3 ) )
% %     
    ax = gca();
    set(ax,axCircle{:})
end

count = 2;
for j = fliplr([14:16,1:4])
    subplot(pY, pX, subNr(count,pX)); hold on
    plot( tNon, sin(2*pi*tNon),':','Color', [1,1,1]*0.1)
    count = count+1;
    
    plot( tNon, -sim.Om0.pFire(:, j,1) ,'k','LineWidth',lw)
    plot(  tNon, sim.Om10.pFire(:, j,1) ,'Color',col(2,:) ,'LineWidth',lw)
%     axis off
    lim = max( sim.Om10.pFire(:, j ,1)); 
    locs = tNon(sim.Om10.spikeInds{j});
    plot( [1,1]'*locs   ,  [lim, 0]'*0.6,'k','LineWidth',slw)
    try
        locs = tNon(sim.Om0.spikeInds{j});
        plot( [1,1]'*locs   ,  [0, -lim]'*0.6 ,'k','LineWidth',slw)
    end
    plot( [-0.5,2],[0,0],'Color',[1,1,1]*0.5)
%     text(0,0,num2str( dT(j)*1e3 ) )
%     
    ax = gca();
    set(ax,axCircle{:})
end



for j = [5,13]
    subplot(pY, pX, reshape( subNr(j-4,2:pX-1),[],1 )); hold on
    plot( tNon, sin(2*pi*tNon),':','Color', [1,1,1]*0.1)
    
    plot( tNon, -sim.Om0.pFire(:, j,1) ,'k','LineWidth',lw)
    plot(  tNon, sim.Om10.pFire(:, j,1) ,'Color',col(2,:) ,'LineWidth',lw)
%     axis off
    lim = max( sim.Om10.pFire(:, j ,1)); 
    locs = tNon(sim.Om10.spikeInds{j});
    plot( [1,1]'*locs   ,  [lim, 0]'*0.6,'k','LineWidth',slw)
    try
        locs = tNon(sim.Om0.spikeInds{j});
        plot( [1,1]'*locs   ,  [0, -lim]'*0.6 ,'k','LineWidth',slw)
    end
    plot( [-0.5,2],[0,0],'Color',[1,1,1]*0.5)
%     text(0.5,0,num2str( dT(j)*1e3 ) )
%     
    ax = gca();
    set(ax,axCircle{:})
end



axCirc = {...
    'PlotBoxAspectRatio',[3 4 4],...
    'DataAspectRatio',[1 1 1],...
    'XLim',[-1,1],...
    'YLim',[-1,1],...
    };
subplot(pY, pX, reshape( subNr(3:7,2:pX-1),[],1) ); hold on
    x = sin(  linspace(0,1,100)*2*pi );
    y = cos( linspace(0,1,100)*2*pi );
    plot(x,y,'k','LineWidth',2)
    axis square


    x = sin(  linspace(0,1,17)*2*pi );
    y = cos( linspace(0,1,17)*2*pi );
%     a = scatter(x,y,50,'ko','filled');

    ax = gca();
    set(ax, axCirc{:})
    axis off

%     for j = 1:16
%         text(  -xyz(2,circlePoints(j))/150*1.3, xyz(3,circlePoints(j))/150*1.3, num2str(  dT(j)*1e3  )  )
%     end
    
% axExplain = {'Visible','off',...
%             'XTick',1:0.5:2,...
%             'XTickLabel',0:0.5:1,...
%             'YTick',[-1,1],...
%             };
% subplot(pY, pX, subNr(1,2:pX-1) )
%     text(1,-1,'Stroke phase (t/T)')
%     text(0.5,1,'Prob. Fire')
%     ax = gca();
%     set(ax,axExplain{:})
   
%% 
figure()
x = sin(  linspace(0,1,100)*2*pi );
y = cos( linspace(0,1,100)*2*pi );
plot(x,y,'k','LineWidth',2)



% angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
% [V,I_sort] = sort(angleDeg,'ascend');
% 
% % Ind = circleIndices(I_sort);
% circlePoints = circleIndices(I_sort); 
% sidePoints = circlePoints( find( mod(V,90) == 0));
dTcor = dT*1e3;
noSpike = find( isnan(dT) )
dTcor( noSpike ) = 10;

for j = 1:16
    yP = xyz(2, circlePoints( j  ) )/150;
    zP = xyz(3, circlePoints( j ) )/150;
    r = yP^2 + zP^2;
    theta = atan2(zP,yP);
    if theta< 0
        theta = theta + 2*pi;
    end
    theta*180/pi;
    
    thetaVec(j) = theta;
    rNew(j) = (r+dTcor(j)/10);
    yNew(j) = cos(theta)*rNew(j);
    zNew(j) = sin(theta)*rNew(j);
    hold on
    scatter( yNew(j), zNew(j), 10,'filled','k')
    
end
rNew = [rNew, rNew(1)];
thetaVec = [thetaVec, thetaVec(1)];
thetaVec(9:end) = thetaVec(9:end) + 2*pi;
yNew = [yNew, yNew(1)];
zNew = [zNew, zNew(1)];
plot(yNew,zNew)

% 

figure()
plot(thetaVec, rNew)

thetaVecInterp = linspace( thetaVec(1), thetaVec(end), 161);

yy = spline(thetaVec,rNew,thetaVecInterp);
plot(thetaVec, rNew,'o',thetaVecInterp,yy)

    
%% Setting paper size for saving 
% 
% set(0,'DefaultAxesFontSize',8)% .
% set(0,'DefaultLegendFontSize',8)% .
% 
% 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% tightfig;
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
