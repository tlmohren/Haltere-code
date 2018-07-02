clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

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

%% apply neural encoding

STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );

tL = size(FEA(1).strain,2);
fSamp = 1000;
subSamp =10;
STAt = linspace(-39,0,40*subSamp);

parse_extra = 5;

calib_param_max = [0.00298239604088481,0.0212443192274280,0.0372780009360082,0.0479350861876651,0.0514646550697229,0.0475135494972593,0.0363905581996612,0.0198191699014800,0.00299958384555981,0.0201797042231435,0.0370669398494605,0.0483996355340012,0.0524251480251299,0.0487363345563017,0.0377914840661745,0.0213792050886826];
% 
for j = 1:length(FEA)
    for k = 1:length(FEA(j).circleInds)
       strainTemp = FEA(j).strain( FEA(j).circleInds(k),:);
        t_plot = (1:tL)/fSamp;
        tNew = linspace(t_plot(1),t_plot(end), tL*subSamp ) ; 
        strainInterp(:,j) =  interp1(t_plot,strainTemp ,tNew,'spline');
       
        STA = STAfun(STAt);
        strainConv = conv( [zeros(1,length(STA)-1),strainInterp(:,j)'], fliplr( STA), 'valid');
        calib_param(j,k) = max(  strainConv );
%         FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param(j,k) );
        FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param_max(k) );
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
selected_dots = 5:9;
len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

axCircle= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)]}; 
axZero= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''}}; 

fig2 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)-[width,0]*100,  width*100, height*100]); %<- Set size

for kl =1:length(selected_dots)
    k = selected_dots(kl);
    subplot(5,1,kl)
    hold on
    for j = 1:length(FEA)
        plot( t_plot,FEA(j).strain(FEA(j).circleInds(k),It )  )
    end
    ylabel('$\epsilon$')
    ax = gca();
    if kl == 5
        set(ax,axCircle{:})
        xlabel('Time (s)')
    else
        set(ax,axZero{:})
    end
end
%% 
fig3 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig3, 'Position', [fig3.Position(1:2)  width*100, height*100]); %<- Set size
subplot(211)
    tSTA = -40:0.1:0;
    plot(tSTA,STAfun(tSTA))
subplot(212)
    s = -1:0.01:1;
    plot(s,NLDfun(s))
    ylabel('blabla')
    %%   

len = 1001;
start = 351;
It10 = start:(start+len-1);
t_plot = (0:len-1)*0.0001;

fig4 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig4, 'Position', [fig4.Position(1:2)+[width,0]*100 width*100, height*100]); %<- Set size
    
for kl =1:length(selected_dots)
    k = selected_dots(kl);
    subplot(5,1,kl)
    hold on
   for j = 1:length(FEA)
        plot( t_plot,FEA(j).pFire( k , It10) )
   end
    ax = gca();
    if kl == 5
        set(ax,axCircle{:})
   xlabel('Time (s)')
    else
        set(ax,axZero{:})
    end
   ylabel('Prob. of Firing')
end
        
%% 

fig5 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig5, 'Position', [fig5.Position(1:2)+[width*2,0]*100 width*100, height*100]); %<- Set size

for kl =1:length(selected_dots)
    k = selected_dots(kl);
    subplot(5,1,kl)
   for j = 1:length(FEA)
    hold on
        spike_ind_temp = FEA(j).spikeInds{k};
        which_ = spike_ind_temp<=It10(end);
            try
               plot( [1,1]'*t_plot(  spike_ind_temp(which_) - start ),...
                       [8-j,7.2-j]'* ones(1, length(  spike_ind_temp(which_) )) ,'k')
              FEA(j).timing{k} = t_plot(  spike_ind_temp(which_) - start );
               if mod(j,2) ==0
                   FEA(j).dt{k} = FEA(j).timing{k} - FEA(j-1).timing{k};
                   text(0.01, 8-j,['$\Delta$ t=',num2str(FEA(j).dt{k}(1))])
               end
            end
   end
    ax = gca();
    if kl == 5
        set(ax,axCircle{:})
    xlabel('Time(s)')
    else
        set(ax,axZero{:})
    end
end

%% 
fig6= figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig6, 'Position', [fig6.Position(1:2)+[width*3,0]*100 width*100, height*100]); %<- Set size
    
for kl =1:length(selected_dots)
    k = selected_dots(kl);
    subplot(5,1,kl)
    hold on
   for j = 1:length(FEA)
        spike_ind_temp = FEA(j).spikeInds{k};
       plot( t_plot,FEA(j).pFire( k , It10) )
       try
           plot( [1,1]'*t_plot(  spike_ind_temp(which_) - start ),...
               [8-j,7.2-j]'* ones(1, length(  spike_ind_temp(which_) ))./8 ,'k')
           FEA(j).timing{k} = t_plot(  spike_ind_temp(which_) - start );
           if mod(j,2) ==0
               FEA(j).dt{k} = FEA(j).timing{k} - FEA(j-1).timing{k};
               text(0.01, (8-j)/8,['$\Delta$ t=',num2str(FEA(j).dt{k}(1))])
           end
       end
   end
%     
    ax = gca();
    if kl == 5
        set(ax,axCircle{:})
   xlabel('Time (s)')
    else
        set(ax,axZero{:})
    end
   ylabel('Prob. of Firing')
end

%% Setting paper size for saving 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure4_spikeStrain'], '-dsvg');
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure4_spikeStrain' ], '-dpng', '-r600');


set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig3,'InvertHardcopy','on');
set(fig3,'PaperUnits', 'inches');
papersize = get(fig3, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig3, 'PaperPosition', myfiguresize);
print(fig3, ['figs' filesep 'Figure3_spikeEncoder' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig3, 'PaperPosition', myfiguresize);
print(fig3, ['figs' filesep 'Figure3_spikeEncoder'], '-dsvg', '-r600');


set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig4,'InvertHardcopy','on');
set(fig4,'PaperUnits', 'inches');
papersize = get(fig4, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig4, 'PaperPosition', myfiguresize);
print(fig4, ['figs' filesep 'Figure4_spikeFire' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig4, 'PaperPosition', myfiguresize);
print(fig4, ['figs' filesep 'Figure4_spikeFire'], '-dsvg');


set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig5,'InvertHardcopy','on');
set(fig5,'PaperUnits', 'inches');
papersize = get(fig5, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig5, 'PaperPosition', myfiguresize);
print(fig5, ['figs' filesep 'Figure4_spikeSpikes'], '-dsvg');
myfiguresize = [left, bottom, width, height];
set(fig5, 'PaperPosition', myfiguresize);
print(fig5, ['figs' filesep 'Figure4_spikeSpikes' ], '-dpng', '-r600');


set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig6,'InvertHardcopy','on');
set(fig6,'PaperUnits', 'inches');
papersize = get(fig6, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig6, 'PaperPosition', myfiguresize);
print(fig6, ['figs' filesep 'Figure4_spikeFireSpikes'], '-dsvg');
myfiguresize = [left, bottom, width, height];
set(fig6, 'PaperPosition', myfiguresize);
print(fig6, ['figs' filesep 'Figure4_spikeFireSpikes' ], '-dpng', '-r600');

