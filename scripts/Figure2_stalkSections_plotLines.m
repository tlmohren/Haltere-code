clc;clear all;close all

addpathFolderStructureHaltere()
run('config_file.m')

loadName = 'FEA_processed_data';
load(['data' filesep loadName],'FEA')

%% deformation in angles 

fig1 = figure();
    width = 3.5;     % Width in inches,   find column width in paper 
    height = 2.5;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2) width*100, height*100]); %<- Set size
    lineType = {'-',':','--'};
    
for j = [1,3] 
    subplot(331); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) ,lineType{j} );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(334 ); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) ,lineType{j}   );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(337); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It)  ,lineType{j} );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end


%%  
for j = 2% 
    subplot(3,3,2 ); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(3,3,5 ); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(3,3,8 ); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end

%% Figure 2

axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15] ,...
               'YLim',[-1,1]*1e-4}; 
           
for j = 4 
    subplot(333); hold on 
        p1 = plot( t_plot, FEA(j).yAngle(It) );
        ylabel( deformLabels{1} );
        ax = gca();
        set(ax,axOpts_dphi{:})
        
    subplot(336); hold on 
        p1 = plot( t_plot, FEA(j).zAngle(It) );
        ylabel( deformLabels{2} );
        ax = gca();
        set(ax,axOpts_dtheta{:})
        
    subplot(339); hold on 
        p1 = plot( t_plot, FEA(j).twistAngle(It) );
        xlabel('Time (s)'); ylabel( deformLabels{3} );
        ax = gca();
        set(ax,axOpts_dgamma{:})
end 
 
%% Setting paper size for saving 
set(gca, 'LooseInset', get(gca(), 'TightInset')); % remove whitespace around figure
tightfig;
set(fig1,'InvertHardcopy','on');
set(fig1,'PaperUnits', 'inches');
papersize = get(fig1, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_deformStalk' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure2_deformStalk'], '-dsvg', '-r600');
