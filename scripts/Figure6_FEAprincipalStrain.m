clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

%%
loadName = 'figure6_strainData';
saveName = 'figure6_strainData';

renew_data_load = false
% renew_data_load = true
outer_dimension = [150,150,243.435,243.435];
if renew_data_load
    FEA(1).name = 'Haltere_CraneFly_Sphere_Om0';
    FEA(2).name = 'Haltere_CraneFly_Sphere_Om10';
    FEA(3).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om0';
    FEA(4).name = 'Haltere_CraneFly_ellipsoidVerCrossStalk_Om10';
    for j =  1:length(FEA)
        tic
        [FEA(j).xyz, FEA(j).strain, ~] = loadCSV( ['data' filesep  FEA(j).name], { 'eXX','eYY','eZZ','eXY','eXZ','eYZ' });        toc 
    end
    circleDistance = 300;               % distance from base to haltere 
    er = 0.01;
    for j = 1:length(FEA)
        sideL = outer_dimension(j); 
        xMatch = find(   abs(abs( FEA(j).xyz(:,1) ) - circleDistance)  <er );
        yMatch = find(   abs(abs( FEA(j).xyz(:,2) ) - sideL) <er & ...
                     abs( FEA(j).xyz(:,3) )   <er);
        zMatch = find(   abs(abs( FEA(j).xyz(:,3) ) - sideL)  <er & ...
                     abs( FEA(j).xyz(:,2) )   <er);
        FEA(j).sideInds = intersect(xMatch,yMatch);
        FEA(j).topInds = intersect(xMatch,zMatch);
    end
    save(['data' filesep saveName],'FEA')
else
    load(['data' filesep loadName],'FEA')
end

%% find principle strains

for j = 1:length(FEA)
    for k = 1:2
       exTop = FEA(j).strain( FEA(j).topInds(k), : , [1] );
       eyTop = FEA(j).strain( FEA(j).topInds(k), : , [2] );
       exyTop = FEA(j).strain( FEA(j).topInds(k), : , [4] );
        
       FEA(j).epTop(k,:) = 0.5*(exTop+eyTop) + sqrt( (0.5*(exTop-eyTop)).^2 + (exyTop).^2 );
       FEA(j).thetapTop(k,:) = 0.5*atan( exyTop./(exTop-eyTop) );
       
       
       exSide = FEA(j).strain( FEA(j).sideInds(k), : , [1] );
       eySide = FEA(j).strain( FEA(j).sideInds(k), : , [3] );
       exySide = FEA(j).strain( FEA(j).sideInds(k), : , [5] );
        
       FEA(j).epSide(k,:) = 0.5*(exSide+eySide) + sqrt( (0.5*(exSide-eySide)).^2 + (exySide).^2 );
       FEA(j).epSideNoShear(k,:) = 0.5*(exSide+eySide) + sqrt( (0.5*(exSide-eySide)).^2 + 0*(exySide).^2 );
       FEA(j).thetapSide(k,:) = 0.5*atan( exySide./(exSide-eySide) );
       FEA(j).thetapSideNoShear(k,:) = 0.5*atan( FEA(j).strain( FEA(3).sideInds(k), : , [5] )./(exSide-eySide) );
    end
end

%% order of magnitude comparison 

for j = [2,4]
    figure()
    for k = 1:2
        subplot(411); hold on 
        plot( FEA(j).strain(FEA(j).sideInds(k), : , [1] ) ) 
        subplot(412); hold on 
        plot( FEA(j).strain(FEA(j).sideInds(k), : , [3] )  )
        subplot(413); hold on 
        plot( FEA(j).strain(FEA(j).sideInds(k), : , [5] ) )


        subplot(414); hold on 
        plot(  FEA(j).epSide(k,:) )

    end
end

for j = [1,3]
    figure()
    for k = 1:2
        subplot(411); hold on 
        plot( FEA(j).strain(FEA(j).sideInds(k), : , [1] ) ) 
        subplot(412); hold on 
        plot( FEA(j).strain(FEA(j).sideInds(k), : , [3] )  )
        subplot(413); hold on 
        plot( FEA(j).strain(FEA(j).sideInds(k), : , [5] ) )


        subplot(414); hold on 
        plot(  FEA(j).epSide(k,:) )

    end
end

%% 
figure()
    for j = 1:2
%     for k = 1:2
        n_plot = 2*(j-1)+k
        subplot(2,2,n_plot)
        plot( squeeze(  FEA(j+2).strain(FEA(j+2).sideInds(1), : ,: ) )  )
%     end
        legend({ 'eXX','eYY','eZZ','eXY','eXZ','eYZ' })
        axis([0,200,-2e-4,2e-4])
    end
  
%% 
figure()
    subplot(211)
    plot( FEA(3).strain( FEA(3).sideInds(1),:,5) )
    hold on
    plot( FEA(4).strain( FEA(4).sideInds(1),:,5) )
%     legend('$\tau_{xz, \text{flapping}}$', '$\tau_{xz, \text{flap \& rot}}$','Interpreter','latex')
    legend('$\tau_{xz~~flapping}$','$\tau_{xz~~flap~and~rot}$','location','NorthWest')
    xlabel('Time index'); ylabel('$\tau_{xz}$')
    title('shear strain at side of halteres, with cross-shaped stalk')
     subplot(212) 
    plot( FEA(3).strain( FEA(3).sideInds(1),:,5) - FEA(4).strain( FEA(4).sideInds(1),:,5) )
    legend('$\tau_{xz~~flapping} - \tau_{xz~~flap~and~rot}$','location','NorthWest')
    hold on
    xlabel('Time index'); ylabel('$\tau_{xz}$')
%     for j = 1:2
% %     for k = 1:2
%         n_plot = 2*(j-1)+k
%         subplot(2,2,n_plot)
%         plot( squeeze(  FEA(j+2).strain(FEA(j+2).sideInds(1), : ,: ) )  )
% %     end
%         legend({ 'eXX','eYY','eZZ','eXY','eXZ','eYZ' })
%         axis([0,200,-2e-4,2e-4])
%     end
        
%     
%     
% %% 
%     figure()
% for j = [4]
%     for k = 1:2
%         subplot(211); hold on 
%         title('principal strain without shear')
%         plot(  FEA(j).epSide(k,:) )
%         subplot(212); hold on 
%         plot(  FEA(j).epSideNoShear(k,:) )
%         title('principal strain with shear')
% 
%     end
% end
% 
%     figure()
% for j = [4]
%     for k = 1:2
%         subplot(211); hold on 
%         title('principal strain without shear')
%         plot(  FEA(j).thetapSide(k,:) )
%         subplot(212); hold on 
%         plot(  FEA(j).thetapSideNoShear(k,:) )
%         title('principal strain with shear')
% 
%     end
% end
% 
% 
% %% 
% axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''}}; 
% axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] }; 
% 
% len = 101;
% start = 35;
% It = start:(start+len-1);
% t_plot = (0:len-1)*0.001;
% 
% 
%     for j = 1:length(FEA)
% fig1 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig1, 'Position', [fig1.Position(1:2)-[width,0]*100 width*100, height*100]); %<- Set size
%         subplot(211);hold on 
%              plot(t_plot, FEA(j).epTop(:,It)')
%             ax = gca();
%             set(ax,axOptsStrainTop{:})
%         subplot(212); hold on 
%             plot(t_plot, FEA(j).thetapTop(:,It)')
%             ax = gca();
%             set(ax,axOptsStrainTop{:})
%     end
%     for j = 1:length(FEA)
% fig2 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size
%         subplot(211); hold on 
%             plot(t_plot, FEA(j).epSide(:,It)')
%             ax = gca();
%             set(ax,axOptsStrainSide{:})
% 
%         subplot(212); hold on 
%             plot(t_plot, FEA(j).thetapSide(:,It)')
%             ax = gca();
%             set(ax,axOptsStrainSide{:})
%     end
% 
% %% deformation in angles 
% 
% fig1 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
% len = 101;
% start = 35;
% It = start:(start+len-1);
% t_plot = (0:len-1)*0.001;
% 
% axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2.5e-3,'YTick',[-1,0,1]*2e-3}; 
% axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 
% 
% for j = [1,3]
%     subplot(211); hold on 
%         plot(t_plot, FEA(j).strain( FEA(j).topInds, It) )
%         ax = gca();
%         set(ax,axOptsStrainTop{:})
%         ylabel('$\epsilon$','rotation',0)
%     subplot(212); hold on 
%         plot(t_plot, FEA(j).strain( FEA(j).sideInds, It) )
%         ax = gca();
%         set(ax,axOptsStrainSide{:})
%         xlabel('Time (s)')
%         ylabel('$\epsilon$','rotation',0)
% end
% 
% %% Setting paper size for saving 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig1,'InvertHardcopy','on');
% set(fig1,'PaperUnits', 'inches');
% papersize = get(fig1, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig1, 'PaperPosition', myfiguresize);
% print(fig1, ['figs' filesep 'Figure6_strainPlotOM0' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig1, 'PaperPosition', myfiguresize);
% print(fig1, ['figs' filesep 'Figure6_strainPlotOM0'], '-dsvg', '-r600');
% % 
% 
% 
% 
% 
% %% deformation in angles 
% 
% fig2 = figure();
%     width = 2;     % Width in inches,   find column width in paper 
%     height = 2;    % Height in inches
%     set(fig2, 'Position', [fig2.Position(1:2) width*100, height*100]); %<- Set size
% len = 101;
% start = 35;
% It = start:(start+len-1);
% t_plot = (0:len-1)*0.001;
% 
% axOptsStrainTop = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)],'XTickLabel',{'','',''} ,'YLim',[-1,1]*2.5e-3,'YTick',[-1,0,1]*2e-3}; 
% axOptsStrainSide = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,'YLim',[-1,1]*1.5e-4}; 
% 
% for j = [2,4]
%     subplot(211); hold on 
%         plot(t_plot, FEA(j).strain( FEA(j).topInds, It) )
%         ax = gca();
%         set(ax,axOptsStrainTop{:})
%         ylabel('$\epsilon$','rotation',0)
%     subplot(212); hold on 
%         plot(t_plot, FEA(j).strain( FEA(j).sideInds, It) )
%         ax = gca();
%         set(ax,axOptsStrainSide{:})
%         xlabel('Time (s)')
%         ylabel('$\epsilon$','rotation',0)
% end
% 
% %% Setting paper size for saving 
% set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
% set(fig2,'InvertHardcopy','on');
% set(fig2,'PaperUnits', 'inches');
% papersize = get(fig2, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(fig2, 'PaperPosition', myfiguresize);
% print(fig2, ['figs' filesep 'Figure6_strainPlotOM10' ], '-dpng', '-r600');
% stupid_ratio = 15/16;
% myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
% set(fig2, 'PaperPosition', myfiguresize);
% print(fig2, ['figs' filesep 'Figure6_strainPlotOM10'], '-dsvg', '-r600');
% % 