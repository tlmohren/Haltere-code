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
angleDeg = rad2deg(angle)-90;
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

fsz = 12;
           
set(0,'DefaultAxesFontSize',fsz)% .
set(0,'DefaultLegendFontSize',fsz)% .

%% general figure settings 
[tI,nS,s ] = size( sim.Om0.strain);
t = ( (1:tI) - tI+50) /25;

axOpts = {'XLim',[-0.5,2],...
    'XTick',[0,1,2],...
    'Xgrid','on',...
    'Ygrid','on',...
    'XTickLabel',{'','',''},...
    };
axOptsBottom = {'XLim',[-0.5,2],...
    'XTick',[0,1,2],...
    'Xgrid','on',...
    'Ygrid','on',...
    };

YaxSides = {'YLim',[-1,1]*2e-3     };
YaxTop = {'YLim',[-1,1]*1e-4     };

col = linspecer(2);



%% 
fig1 = figure();
set(fig1, 'Position', [100 300 width*100, height*100]); %<- Set size

pX = 6;
pY = 4;
subGrid = reshape(1:(pX*pY),pX,pY)';
order = [1,3,2,4]; 
subInd = Ind( 1:4:13 );
ySel = [1,1,2,2];
dotCol = {'ob','ob','or','or'}; 

tLW = 2;
ColorlineOpts = {'-.','Color', col(1,:),'LineWidth',tLW; ....
            '-.','Color', col(1,:),'LineWidth',tLW; ....
            '-.','Color', col(2,:),'LineWidth',tLW; ....
            '-.','Color', col(2,:),'LineWidth',tLW; ....
            };

for j = 1:4 
    subplot( pY, pX, subGrid(j,3:end) ); hold on
    plot( t, sim.Om0.strain(:, sidePoints( order(j) ),1) , ColorlineOpts{j,:})
    ax = gca();
    if j ==4
        set(ax,axOptsBottom{:})
    else
        set(ax,axOpts{:})
    end
    if any( order(j) == [1,3] )
        set(ax,YaxSides{:})
    else
        set(ax,YaxTop{:})
    end
end
    

axCirc = {...
    'PlotBoxAspectRatio',[3 4 4],...
    'DataAspectRatio',[1 1 1],...
    'XLim',[-1.2,1.2],...
    'YLim',[-1.2,1.2],...
    };
for j = 1:4
   subplot( pY, pX, subGrid(j,1) ) 
    x = sin(  linspace(0,1,100)*2*pi );
    y = cos( linspace(0,1,100)*2*pi );
    plot(x,y,'k')
    hold on
    
%     axis square
%     axis off
    
    yP = xyz(2, subInd( order(j)  ) )/150;
    zP = xyz(3, subInd( order(j) ) )/150;
    plot( [0,0],[0,1.2],'k')
    plot( [-1.2,0],[0,0],'k')
    text( -1.7,0,'x')
    text( -0.2,1.6, 'z')
    
    scatter( yP, zP, dotCol{j},'filled')
    
    ax = gca();
    set(ax, axCirc{:})
    axis off
end
    
    
    %% 
fig2 = figure();
set(fig2, 'Position', [500 300 width*100, height*100/2]); %<- Set size

pX = 6;
pY = 2;
subGrid = reshape(1:(pX*pY),pX,pY)';

BlacklineOpts = {'-.','Color', 'k','LineWidth',tLW; ....
            '-.','Color', 'k','LineWidth',tLW; ....
            '-.','Color','k','LineWidth',tLW; ....
            '-.','Color', 'k','LineWidth',tLW; ....
            };

for j = 1:4 
    subplot( pY, pX, subGrid( ySel(j) ,3:end) ); hold on
    plot( t, sim.Om0.strain(:, sidePoints( order(j) ),1) , BlacklineOpts{j,:})
    ax = gca();
    if any(j == [3,4] ) 
        set(ax,axOptsBottom{:})
        
    else
        set(ax,axOpts{:})
    end
    if any( order(j) == [1,3] )
        set(ax,YaxSides{:})
    else
        set(ax,YaxTop{:})
    end
end
    
for j = 1:4
   subplot( pY, pX, subGrid(ySel(j),1) ) 
    x = sin(  linspace(0,1,100)*2*pi );
    y = cos( linspace(0,1,100)*2*pi );
    plot(x,y,'k')
    hold on
    
    axis square
    axis off
    
    yP = xyz(2, subInd( order(j)  ) )/150;
    zP = xyz(3, subInd( order(j) ) )/150;
    scatter( yP, zP, dotCol{j},'filled')
end
    
    







%% 
fig3 = figure();
set(fig3, 'Position', [ 900, 300, width*100, height*100]); %<- Set size

mLW = 1.5;
ColorlineOpts = {'-','Color', col(1,:),'LineWidth',mLW; ....
            '-','Color', col(1,:),'LineWidth',mLW; ....
            '-','Color', col(2,:),'LineWidth',mLW; ....
            '-','Color', col(2,:),'LineWidth',mLW; ....
            };

        
        
pX = 6;
pY = 4;
subGrid = reshape(1:(pX*pY),pX,pY)';

for j = 1:4 
    subplot( pY, pX, subGrid(j,3:end) ); hold on
    plot( t, sim.Om10.strain(:, sidePoints( order(j) ),1) , ColorlineOpts{j,:})
    ax = gca();
    if j ==4
        set(ax,axOptsBottom{:})
    else
        set(ax,axOpts{:})
    end
    if any( order(j) == [1,3] )
        set(ax,YaxSides{:})
    else
        set(ax,YaxTop{:})
    end
end
    
for j = 1:4
   subplot( pY, pX, subGrid(j,1) ) 
    x = sin(  linspace(0,1,100)*2*pi );
    y = cos( linspace(0,1,100)*2*pi );
    plot(x,y,'k')
    hold on
    
    axis square
    axis off
    
    yP = xyz(2, subInd( order(j)  ) )/150;
    zP = xyz(3, subInd( order(j) ) )/150;
    scatter( yP, zP, dotCol{j},'filled')
end 



    
    %% 
fig4 = figure();
set(fig4, 'Position', [1300 300 width*100, height*100/2]); %<- Set size

pX = 6;
pY = 2;
subGrid = reshape(1:(pX*pY),pX,pY)';


pX = 6;
pY = 2;
subGrid = reshape(1:(pX*pY),pX,pY)';



for j = 1:4 
    subplot( pY, pX, subGrid( ySel(j) ,3:end) ); hold on
    plot( t, sim.Om0.strain(:, sidePoints( order(j) ),1) , BlacklineOpts{j,:})
    ax = gca();
    if any( j == [3,4])
        set(ax,axOptsBottom{:})
%         j
    else
        set(ax,axOpts{:})
    end
    if any( order(j) == [1,3] )
        set(ax,YaxSides{:})
    else
        set(ax,YaxTop{:})
    end
end
    
for j = 1:4
   subplot( pY, pX, subGrid(ySel(j),1) ) 
    x = sin(  linspace(0,1,100)*2*pi );
    y = cos( linspace(0,1,100)*2*pi );
    plot(x,y,'k')
    hold on
    
    axis square
    axis off
    
    yP = xyz(2, subInd( order(j)  ) )/150;
    zP = xyz(3, subInd( order(j) ) )/150;
    scatter( yP, zP, dotCol{j},'filled')
end
    
for j = 1:4 
    subplot( pY, pX, subGrid( ySel(j) ,3:end) ); hold on
    plot( t, sim.Om10.strain(:, sidePoints( order(j) ),1) , ColorlineOpts{j,:})
    ax = gca();
end
    




    
%% Setting paper size for saving 

% % % % set(0,'DefaultAxesFontSize',8)% .
% % % % set(0,'DefaultLegendFontSize',8)% .
% % % % 
% % % % 
% % % % set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% % % % tightfig;
% % % % % % Here we preserve the size of the image when we save it.
% % % % 
% % % % 
% % % % set(fig1,'InvertHardcopy','on');
% % % % set(fig1,'PaperUnits', 'inches');
% % % % papersize = get(fig1, 'PaperSize');
% % % % left = (papersize(1)- width)/2;
% % % % bottom = (papersize(2)- height)/2;
% % % % myfiguresize = [left, bottom, width, height];
% % % % set(fig1, 'PaperPosition', myfiguresize);
% % % % 
% % % % % Saving figure 
% % % % print(fig1, ['figs' filesep 'Figure_StrainsFR' ], '-dpng', '-r600');
% % % % 
% % % % % total hack, why does saving to svg scale image up???
% % % % stupid_ratio = 15/16;
% % % % myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% % % % set(fig1, 'PaperPosition', myfiguresize);
% % % % 
% % % % print(fig1, ['figs' filesep 'Figure_StrainsFR' ], '-dsvg');
