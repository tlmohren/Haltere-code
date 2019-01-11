clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

loadName = 'FEA_processed_data';
load(['data' filesep loadName],'FEA')

%% apply neural encoding

FEA(13).circleInds = FEA(1).circleInds;
FEA(14).circleInds = FEA(2).circleInds;

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

selected_dots = 5:9;
len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

axCircle= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)]}; 
axZero= { 'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''}}; 
axCircleSpike= { 'YLim',[0,1.6],'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)]}; 
axZeroSpike= { 'YLim',[0,1.6],'XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''}}; 

calib_param_cross = [0.0047, 0.0803, 0.0047, 0.0803]; 
for jj = 1:length(FEA)/2  
%     for j = 1:2%length(FEA)
    for j = [1,2]+(2*(jj-1)) 
        for k = 1:length(FEA(j).circleInds)
         tL = size(FEA(j).strain,2);
            strainTemp = FEA(j).strain( FEA(j).circleInds(k),:);
            t_plot = (1:tL)/fSamp;
            tNew = linspace(t_plot(1),t_plot(end), tL*subSamp ) ; 
            strainInterp =  interp1(t_plot,strainTemp ,tNew,'spline');

            STA = STAfun(STAt);
            strainConv = conv( [zeros(1,length(STA)-1),strainInterp], fliplr( STA), 'valid');
%             strainConv = conv( [zeros(1,length(STA)-1),strainInterp], ( STA), 'valid');
            calib_param(j,k) = max(  strainConv );
    %         FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param(j,k) );
%             if jj ~=7
%                 j
                FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param_max(k) );
                FEA(j).spikeInds{k} = findSpikes( FEA(j).pFire(k,:) ); 
%             elseif jj == 7
%                 
%                 FEA(j).pFire(k,:) = NLDfun( strainConv/calib_param_cross(k) );
%                 FEA(j).spikeInds{k} = findSpikes( FEA(j).pFire(k,:) ); 
%             end
        end
    end
    for k = 1:length(FEA(j).circleInds)
       FEA(jj*2).dI = FEA(jj*2).spikeInds{k}(end) - FEA(jj*2-1).spikeInds{k}(end);
       FEA(jj*2).dT(k) = abs(FEA(jj*2).dI/1e4 );
    end 
end
 
%% figure prep 

len1 = 26;
start1 = 166;
It1 = start1:(start1+len1-1);
t_plot1 = linspace(0,25,len1);

len10 = 261;
start10 = 1651;
It10 = start10:(start10+len10-1);
t_plot10 = linspace(0,25,len10);
 
spike_order = [5:13, fliplr([13:16, 1:5]) ];
%% 

fig7= figure();
    width = 5;     % Width in inches,   find column width in paper 
    height = 5;    % Height in inches
    set(fig7, 'Position', [fig7.Position(1:2)-[width,height]*100 width*100, height*100]); %<- Set size

subplot(1,2,1); hold on 
    scatter( FEA(13).xyz( FEA(13).circleInds(5:13),2 ),  FEA(13).xyz( FEA(13).circleInds(5:13),3 ) +150+350,'k','filled')
    scatter( FEA(13).xyz( FEA(13).circleInds( [13:16, 1:5] ),2 ),  FEA(13).xyz( FEA(13).circleInds( [13:16, 1:5] ),3 ) +150 ,'k','filled')
    theta = 0:0.001:pi;
    plot( -sin(theta)*150,cos(theta)*150 +150+350,'k')
    plot( sin(theta)*150,cos(theta)*150 +150,'k')
    axis equal
    axis([-180,180,-10,900])
    axis off

subplot(122) ; hold on 
%     rectangle('Position',[1,2.9,24,4.8],'Curvature',0,'FaceColor',[1,1,1]*0.95)
%     rectangle('Position',[1,-2.6,24,4.8],'Curvature',0,'FaceColor',[1,1,1]*0.95)
    plot(  t_plot1, FEA(13).phi(1,It1) +10 ,'k')

for kl = 1:length(spike_order)
    k = spike_order( kl ) ;
    spike_ind_temp = FEA(13).spikeInds{k};
        which_ = (spike_ind_temp<=It10(end)) & (spike_ind_temp>=It10(1));

    if any(spike_ind_temp) & kl <=9
       plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10+1 ),...
           ([9 ; 8.6] - kl*0.5 -1)  * ones(1, length(  spike_ind_temp(which_) )),'k')
    elseif any(spike_ind_temp) 
       plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10+1 ),...
           ([9 ; 8.6] - kl*0.5 - 2)  * ones(1, length(  spike_ind_temp(which_) )),'k')
    end
    
    spike_ind_temp = FEA(14).spikeInds{k};
        which_ = (spike_ind_temp<=It10(end)) & (spike_ind_temp>=It10(1));

    if any(spike_ind_temp) & kl <=9
       plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10+1 ),...
           ([9 ; 8.6] - kl*0.5 -1)  * ones(1, length(  spike_ind_temp(which_) )),'r')
    elseif any(spike_ind_temp) 
       plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 +1),...
           ([9 ; 8.6] - kl*0.5 - 2)  * ones(1, length(  spike_ind_temp(which_) )),'r')
    end
%     end
    if kl == 1
           plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 +1 ),...
               ([-3,10])  * ones(1, length(  spike_ind_temp(which_) )),':k')
    elseif kl == 9
           plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 +1),...
               ([-3,10])'  * ones(1, length(  spike_ind_temp(which_) )),':k')
    end
end
 
xlabel('Time (ms)')
ylabel('Flapping angle $\phi(t)$') 
set(gca(),'YTick', [10-pi/2,10,10+pi/2] ,'YTickLabel',{'$\frac{\pi}{2}$',0,'$\frac{\pi}{2}$' } )
title('Spikes along circumference')


%% save spike figure 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig7,'InvertHardcopy','on');
set(fig7,'PaperUnits', 'inches');
papersize = get(fig7, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig7, 'PaperPosition', myfiguresize);
print(fig7, ['figs' filesep 'Figure5_spikeSphere' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig7, 'PaperPosition', myfiguresize);
print(fig7, ['figs' filesep 'Figure5_spikeSphere'  ], '-dsvg', '-r600');
