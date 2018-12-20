clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

loadName = 'FEA_processed_data';
load(['data' filesep loadName],'FEA')

%% deformation in angles 
fig3 = figure();
    width = 6;     % Width in inches,   find column width in paper 
    height = 3;    % Height in inches
    set(fig3, 'Position', [fig3.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
 

axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*1e-5};   
lineCols = linspecer(3); 
colIt = 1;
for k = [1,3,4]
    for j = [1,2]
        jk = 2*(k-1)+j; 
        subplot(3,5, j ); hold on 
            p1 = plot( t_plot, FEA(jk).yAngle(It)  ,'color',lineCols(colIt,:) );
            ylabel( deformLabels{1} );
            ax = gca();
            set(ax,axOpts_dphi{:})

        subplot(3,5, j+5 ); hold on 
            p1 = plot( t_plot, FEA(jk).zAngle(It) ,'color',lineCols(colIt,:) );
            ylabel( deformLabels{2} );
            ax = gca();
            set(ax,axOpts_dtheta{:})

        subplot(3,5, j+ 10 ); hold on 
            p1 = plot( t_plot, FEA(jk).twistAngle(It)  ,'color',lineCols(colIt,:) );
            xlabel('Time (s)'); ylabel( deformLabels{3} );
            ax = gca();
            set(ax,axOpts_dgamma{:})
    end
        colIt = colIt +1;
end

           
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*5e-4};   

for k = 5:7
    
    for j = [1,2]
        jk = 2*(k-1)+j;
        subplot(3,5, j +2); hold on 
            p1 = plot( t_plot, FEA(jk).yAngle(It)  ,'color',lineCols(k-4,:) );
            ylabel( deformLabels{1} );
%             ax = gca();
%             set(ax,axOpts_dphi{:})

        subplot(3,5, j+7 ); hold on 
            p1 = plot( t_plot, FEA(jk).zAngle(It) ,'color',lineCols(k-4,:) );
            ylabel( deformLabels{2} );
            ax = gca();
            set(ax,axOpts_dtheta{:})

        subplot(3,5, j+12 ); hold on 
            p1 = plot( t_plot, FEA(jk).twistAngle(It)  ,'color',lineCols(k-4,:) );
            xlabel('Time (s)'); ylabel( deformLabels{3} );
            ax = gca();
            set(ax,axOpts_dgamma{:})
    end
end


axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:t_plot(end)] ,...
               'YLim',[-1,1]*3e-5};   

for k = 5:7
    for j = [ 2]
        jk = 2*(k-1)+j;
        subplot(3,5, 15); hold on 
            p1 = plot( t_plot, FEA(jk).twistAngle(It)  ,'color',lineCols(k-4,:) );
            xlabel('Time (s)'); ylabel( deformLabels{3} );
            ax = gca();
            set(ax,axOpts_dgamma{:})
    end
end

%% Setting paper size for saving 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
set(fig3,'InvertHardcopy','on');
set(fig3,'PaperUnits', 'inches');
papersize = get(fig3, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig3, 'PaperPosition', myfiguresize);

print(fig3, ['figs' filesep 'Figure2_deformBulb' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig3, 'PaperPosition', myfiguresize);
print(fig3, ['figs' filesep 'Figure2_deformBulb'], '-dsvg', '-r600');


%% 
figure();
    subplot(311) 
    plot(  FEA(10).yAngle  -FEA(4).yAngle  )
    subplot(312)
    plot(  FEA(12).yAngle  -FEA(6).yAngle  )
    subplot(313)
    plot(  FEA(14).yAngle(1:150)  -FEA(8).yAngle(1:150)   )