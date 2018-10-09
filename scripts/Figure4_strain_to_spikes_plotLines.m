clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure4_strainData';
saveName = 'figure4_strainData';

renew_data_load = false
% renew_data_load = true
if renew_data_load
%     FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
%     FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(1).name = 'Haltere_CraneFlyLowDensity_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFlyLowDensity_Sphere_Om10';
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
        [~, FEA(j).phi, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'flapangle'});
        
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

fSamp = 1000;
subSamp =10;
STAt = linspace(-39,0,40*subSamp);

parse_extra = 5;
% [2.38819645997784e-05,0.000150857690670539,0.000272864784543775,0.000353816252497571,0.000384130393993237,0.000355476077108979,0.000272961111255132,0.000149362928687585,2.40516028295196e-05,0.000149596181836659,0.000269840834412852,0.000352369735988511,0.000382519062467057,0.000353733505263798,0.000271346349615106,0.000148675107229601]
calib_param_max_old = [0.00298239604088481,0.0212443192274280,0.0372780009360082,0.0479350861876651,0.0514646550697229,0.0475135494972593,0.0363905581996612,0.0198191699014800,0.00299958384555981,0.0201797042231435,0.0370669398494605,0.0483996355340012,0.0524251480251299,0.0487363345563017,0.0377914840661745,0.0213792050886826];
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
for jj = 1%1:length(FEA)/2 
    
    
%     jj = 2;
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
    if jj ~= 7
        selected_dots = 5:9;
    else
        selected_dots = 2:3;
    end
    len = 101;
    start = 35;
    It = start:(start+len-1);
    t_plot = (0:len-1)*0.001;
    
    fig2 = figure();
        width = 2;     % Width in inches,   find column width in paper 
        height = 3;    % Height in inches
        set(fig2, 'Position', [fig2.Position(1:2),  width*100, height*100]); %<- Set size

    for kl =1:length(selected_dots)
%         if 
        k = selected_dots(kl);
        subplot(5,1,kl)
        hold on
        for j = [1,2]+(2*(jj-1))
            plot( t_plot,FEA(j).strain(FEA(j).circleInds(k),It )  )
        end
        if kl == 3
            ylabel('$\epsilon$')
        end
        ax = gca();
        if kl == 5
            set(ax,axCircle{:})
            xlabel('Time (s)')
        else
            set(ax,axZero{:})
        end
    end
   
    % 
    set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
    set(fig2,'InvertHardcopy','on');
    set(fig2,'PaperUnits', 'inches');
    papersize = get(fig2, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(fig2, 'PaperPosition', myfiguresize);
    print(fig2, ['figs' filesep 'Figure3_spikeStrain_'  FEA(j).name(18:end-5) ], '-dpng', '-r600');
    stupid_ratio = 15/16;
    myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
    set(fig2, 'PaperPosition', myfiguresize);
    print(fig2, ['figs' filesep 'Figure3_spikeStrain_'  FEA(j).name(18:end-5) ], '-dsvg', '-r600');
    
%% 
    len10 = 1001;
    start10 = 351;
    It10 = start10:(start10+len10-1);
    t_plot10 = (0:len10-1)*0.0001;
    fig6= figure();
        width = 2;     % Width in inches,   find column width in paper 
        height = 3;    % Height in inches
        set(fig6, 'Position', [fig6.Position(1:2)+[width*1,0]*100 width*100, height*100]); %<- Set size

    for kl =1:length(selected_dots)
        k = selected_dots(kl);
        subplot(5,1,kl)
        hold on
        for j = [1,2]+(2*(jj-1))
            spike_ind_temp = FEA(j).spikeInds{k};
                which_ = spike_ind_temp<=It10(end);
           plot( t_plot10,FEA(j).pFire( k , It10) )
%            try
            if any(spike_ind_temp)
               plot( [1,1]'*t_plot10(  spike_ind_temp(which_) - start10 ),...
                   ([1.75 ; 1.55] - (j-(2*(jj-1)))/4 )  * ones(1, length(  spike_ind_temp(which_) )),'k')
               FEA(j).timing{k} = t_plot10(  spike_ind_temp(which_) - start10 );
               try
                   if mod(j,2) == 0
                       FEA(j).dt{k} = FEA(j).timing{k} - FEA(j-1).timing{k};
                       text(0.01, 1,['$\Delta$ t=',num2str(FEA(j).dt{k}(1))])
                   end
               end
            else
                   text(0.01, 1,['$\Delta$ t= NaN'])
            end
       end
    %     
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

    set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
    set(fig6,'InvertHardcopy','on');
    set(fig6,'PaperUnits', 'inches');
    papersize = get(fig6, 'PaperSize');
    left = (papersize(1)- width)/2;
    bottom = (papersize(2)- height)/2;
    myfiguresize = [left, bottom, width, height];
    set(fig6, 'PaperPosition', myfiguresize);
    print(fig6, ['figs' filesep 'Figure3_spikes_'  FEA(j).name(18:end-5) ], '-dpng', '-r600');
    stupid_ratio = 15/16;
    myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
    set(fig6, 'PaperPosition', myfiguresize);
    print(fig6, ['figs' filesep 'Figure3_spikes_'  FEA(j).name(18:end-5) ], '-dsvg', '-r600');

    
end
