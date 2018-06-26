clc;clear all;close all

addpathFolderStructureHaltere()
% run('config_file.m')

%%
loadName = 'figure4_strainData';
saveName = 'figure4_strainData';

renew_data_load = false
% renew_data_load = true
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
        
        
        rMatch = find( sqrt(FEA(j).xyz(:,2).^2 + FEA(j).xyz(:,3).^2)  >= circleRadius*0.99 );

        yMatch = find( round( abs( FEA(j).xyz(:,2) ), 7) == circleRadius );
        zMatch = find( round( abs( FEA(j).xyz(:,3) ), 7) == circleRadius );

        FEA(j).circleIndsUnsorted = intersect( xMatch,rMatch); 
        
        angle = atan2( FEA(j).xyz( FEA(j).circleIndsUnsorted,3), ...
            FEA(j).xyz( FEA(j).circleIndsUnsorted,2) );
        angleDeg = rad2deg(angle);
        angleDeg(angleDeg<0) = angleDeg(angleDeg<0)+360;
        [V,I_sort] = sort(angleDeg,'ascend');
% 
        FEA(j).circleInds= FEA(j).circleIndsUnsorted(I_sort);

        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end


% 
%% apply neural encoding


STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );

% [tL, n, nS] = size( sim.Om0.strain  )
% [tL, n, nS] = size( [ FEA(1).strain(FEA(1).sideInds(),FEA(1).strain(FEA(1).topInds() ] ,:) )
tL = size(FEA(1).strain,2);
fSamp = 1000;
subSamp =10;
STAt = linspace(-39,0,40*subSamp);

parse_extra = 5;

% 
for j = 1:length(FEA)
    for k = 1:length(FEA(j).circleInds)
       strainTemp = FEA(j).strain( FEA(j).circleInds(k),:);
        t = (1:tL)/fSamp;
        tNew = linspace(t(1),t(end), tL*subSamp ) ; 
        strainInterp(:,j) =  interp1(t,strainTemp ,tNew,'spline');
       
        STA = STAfun(STAt);
        strainConv = conv( [zeros(1,length(STA)-1),strainInterp(:,j)'], fliplr( STA), 'valid');
        calib_param(j,k) = max(  strainConv );
        FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param(j,k) );
        
        
        FEA(j).spikeInds{k} = findSpikes( FEA(j).pFire(k,:) ); 
    end
end

for j = 1:length(FEA)/2
    for k = 1:length(FEA(j).circleInds)
       FEA(j*2).dI = FEA(j*2).spikeInds{k}(end) - FEA(j*2-1).spikeInds{k}(end);
       FEA(j*2).dT(k) = abs(FEA(j*2).dI/1e4 );
    end
end



%% make spike figure
% plot strains 

selected_dots = 5:9;

axCircle= {
    'XLim',[0,0.15],...
%     'Visible','off',...
%     'YLim',[-1,1],...
%     'Ygrid','on',...
%     'XTick',0:1:2,...
}; 

axZero= {
    'XLim',[0,0.15],...
    'XTick',[],...
%     'Visible','off',...
%     'YLim',[-1,1],...
%     'Ygrid','on',...
}; 

% 
%     x = sin(  linspace(0,1,100)*2*pi )*150;
%     y = cos( linspace(0,1,100)*2*pi )*150;
% fig1 = figure();
%     width = 3;     % Width in inches,   find column width in paper 
%     height = 3;    % Height in inches
%     set(fig1, 'Position', [fig1.Position(1:2) - [500,500] width*100, height*100]); %<- Set size
%     hold on
%     plot(x,y,'k','LineWidth',2)
% for kl = 1:length(selected_dots)
%     k = selected_dots(kl);
%     scatter( ...
%     FEA(1).xyz( FEA(1).circleInds(k), 2),...
%     FEA(1).xyz( FEA(1).circleInds(k), 3),...
%     50,'filled')
% end

fig2 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2) - [200,500] width*100, height*100]); %<- Set size


for kl =1:length(selected_dots)
    k = selected_dots(kl);
    subplot(5,1,kl)
    hold on
    for j = 1:length(FEA)
        plot( t,FEA(j).strain(FEA(j).circleInds(k),: )  )
    end
    ylabel('blabla')
    ax = gca();
    if kl == 5
        set(ax,axCircle{:})
    else
        set(ax,axZero{:})
    end
end
%% 
fig3 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig3, 'Position', [fig3.Position(1:2) - [-100,500]  width*100, height*100]); %<- Set size
subplot(211)
    tSTA = -40:0.1:0;
    plot(tSTA,STAfun(tSTA))
subplot(212)
    s = -1:0.01:1;
    plot(s,NLDfun(s))

    ylabel('blabla')
    %%   



fig4 = figure();
    width = 3;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig4, 'Position', [fig4.Position(1:2) - [-400,500] width*100, height*100]); %<- Set size
    
for kl =1:length(selected_dots)
    k = selected_dots(kl);
    subplot(5,1,kl)
   for j = 1:length(FEA)
        plot( tNew,FEA(j).pFire( k ,: ) )
    hold on
        scatter( tNew( FEA(j).spikeInds{k} ),...
            ones(size(FEA(j).spikeInds{k}) )*(1-j*0.1) ,'+')
   end
    
    ax = gca();
    if kl == 5
        set(ax,axCircle{:})
    else
%         set(ax,axCircle{:})
        set(ax,axZero{:})
    end
   xlabel('test')
   ylabel('blabla')
end
        



% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig4,'InvertHardcopy','on');
set(fig4,'PaperUnits', 'inches');
papersize = get(fig4, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig4, 'PaperPosition', myfiguresize);
print(fig4, ['figs' filesep 'Figure4_spikeFire'], '-dsvg', '-r600');
% myfiguresize = [left, bottom, width, height];
% set(fig4, 'PaperPosition', myfiguresize);
% print(fig4, ['figs' filesep 'Figure4_spikeFire' ], '-dpng', '-r600');

%% Setting paper size for saving 

% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig1,'InvertHardcopy','on');
% set(fig1,'PaperUnits', 'inches');
% papersize = get(fig1, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig1, 'PaperPosition', myfiguresize);
% print(fig1, ['figs' filesep 'Figure4_spikeCircle' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig1, 'PaperPosition', myfiguresize);
% print(fig1, ['figs' filesep 'Figure4_spikeCircle'], '-dsvg', '-r600');



% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig2,'InvertHardcopy','on');
% set(fig2,'PaperUnits', 'inches');
% papersize = get(fig2, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig2, 'PaperPosition', myfiguresize);
% print(fig2, ['figs' filesep 'Figure4_spikeStrain'], '-dsvg');
% myfiguresize = [left, bottom, width, height];
% set(fig2, 'PaperPosition', myfiguresize);
% print(fig2, ['figs' filesep 'Figure4_spikeStrain' ], '-dpng', '-r600');
% 
% 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig3,'InvertHardcopy','on');
% set(fig3,'PaperUnits', 'inches');
% papersize = get(fig3, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig3, 'PaperPosition', myfiguresize);
% print(fig3, ['figs' filesep 'Figure4_spikeEncoder' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig3, 'PaperPosition', myfiguresize);
% print(fig3, ['figs' filesep 'Figure4_spikeEncoder'], '-dsvg', '-r600');
% 
% 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig4,'InvertHardcopy','on');
% set(fig4,'PaperUnits', 'inches');
% papersize = get(fig4, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig4, 'PaperPosition', myfiguresize);
% print(fig4, ['figs' filesep 'Figure4_spikeFire' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig4, 'PaperPosition', myfiguresize);
% print(fig4, ['figs' filesep 'Figure4_spikeFire'], '-dsvg');