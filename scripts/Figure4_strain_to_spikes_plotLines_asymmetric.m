clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

loadName = 'FEA_processed_data';
load(['data' filesep loadName],'FEA')

FEA(13).circleInds = FEA(1).circleInds;
FEA(14).circleInds = FEA(2).circleInds;
%% Strain plot prep 
selected_dots = 5:9;
len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

axCircle= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'Xgrid','on'}; 
axZero= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''},'Xgrid','on'}; 
 
%% start figure 
fig1 = figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3.5;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2),  width*100, height*100]); %<- Set size


    
    
t_plot = (0:len-1)*0.001;

subplot(6,1,1)
plot(t_plot, FEA(13).phi(It),'k' )
    axPhi= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''},...
        'Xgrid','on','box','off','YTick',[-pi/2,0,pi/2],'YTickLabel',{'$-\pi/2$','0','$\pi/2$'} }; 
    ax = gca();
    set(ax,axPhi{:})
    ylabel('Flapping angle $\phi$ (rad)')

for kl =1:length(selected_dots) 
    k = selected_dots(kl);
    subplot(6,1,kl+1)
    hold on
    for j = [13,14] 
        plot( t_plot,FEA(j).strain(FEA(j).circleInds(k),It )  )
    end
    if kl == 3
        ylabel('Strain')
    end
    ax = gca();
    if kl == 5
        set(ax,axCircle{:})
        xlabel('Time (s)')
    else
        set(ax,axZero{:})
    end
end

%% Save strain figure 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure4_spikeStrain' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure4_spikeStrain' ], '-dsvg', '-r600');



%% apply neural encoding 
len10 = 1001;
start10 = 351;
It10 = start10:(start10+len10-1);
t_plot10 = (0:len10-1)*0.0001;
 
STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );

fSamp = 1000;
subSamp =10;
STAt = linspace(-39,0,40*subSamp);

parse_extra = 5;
 calib_param_max = [0.000705312124793717,0.00521574728317378,0.00916987624184632,0.0117979355154750,0.0126693024102941,0.0116192206073620,0.00884906250717318,0.00484166706412837,0.000716359983875505,0.00495047057398658,0.00903488798016781,0.0118425516863485,0.0128937828713950,0.0119841717461618,0.00928771453961158,0.00524523737029323];

axCircleSpike= { 'YLim',[0,1.6],'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'Xgrid','on'}; 
axZeroSpike= { 'YLim',[0,1.6],'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''},'Xgrid','on'}; 
  
for j = [13,14] 
    for k = 1:length(FEA(j).circleInds)
     tL = size(FEA(j).strain,2);
        strainTemp = FEA(j).strain( FEA(j).circleInds(k),:);
        t_sim = (1:tL)/fSamp;
        tNew = linspace(t_sim(1),t_sim(end), tL*subSamp ) ; 
        strainInterp =  interp1(t_sim,strainTemp ,tNew,'spline');

        STA = STAfun(STAt);
        strainConv = conv( [zeros(1,length(STA)-1),strainInterp], fliplr( STA), 'valid'); 
        calib_param(j,k) = max(  strainConv ); 
            FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param_max(k) );
            FEA(j).spikeInds{k} = findSpikes( FEA(j).pFire(k,:) );  
    end
end 
 
%% plot spike figure
fig2= figure();
    width = 2;     % Width in inches,   find column width in paper 
    height = 3.5;    % Height in inches
    set(fig2, 'Position', [fig2.Position(1:2)+[width*1,0]*100 width*100, height*100]); %<- Set size

subplot(6,1,1)
plot(t_plot, FEA(13).phi(It),'k' ) 
    ax = gca();
    set(ax,axPhi{:})
    ylabel('Flapping angle $\phi$ (rad)')

for kl =1:length(selected_dots)
    k = selected_dots(kl);
    subplot(6,1,kl+1)
    hold on
    for j = [13,14] 
        spike_ind_temp = FEA(j).spikeInds{k};
            which_ = spike_ind_temp<=It10(end);
       plot( t_plot10,FEA(j).pFire( k , It10) )
        if any(spike_ind_temp)
           plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
               ([1.75 ; 1.55] - (j-10)/4 )  * ones(1, length(  spike_ind_temp(which_) )),'k')
           FEA(j).timing{k} = t_plot10(  spike_ind_temp(which_) - start10 );

        end
    end
    ax = gca();
    if kl == 5
        set(ax,axCircleSpike{:})
   xlabel('Time (s)')
    else
        set(ax,axZeroSpike{:})
    end
    if kl == 3
    ylabel('Prob. of Firing')
    end
end

%% save spike figure 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig2,'InvertHardcopy','on');
set(fig2,'PaperUnits', 'inches');
papersize = get(fig2, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure4_spikes '   ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig2, 'PaperPosition', myfiguresize);
print(fig2, ['figs' filesep 'Figure4_spikes '  ], '-dsvg', '-r600');


