%plot neural encodier 
clc;clear all;close all

addpathFolderStructureHaltere()


STAfreq = 0.5;
STAwidth = 5;
STAdelay = 5;
NLDgrad = 20;
NLDshift = 0.8;
[STAfun,NLDfun]=createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift );


 

    STAt = -39:0.1:0;   
    NLDrange = -1:0.01:1;
%     figure(100);
fig1 = figure();
    width = 1.5;     % Width in inches,   find column width in paper 
    height = 2;    % Height in inches
    set(fig1, 'Position', [fig1.Position(1:2)-[width*100,0] width*100, height*100]); %<- Set size
    subplot(211)
        plot(STAt, STAfun(STAt) );
        hold on;
%         drawnow
    xlabel('Time (ms)'); ylabel('Strain')
    subplot(212)
        plot(NLDrange, NLDfun(NLDrange) );
            hold on;
%             drawnow; 
%             grid on

xlabel('$\xi$'); ylabel('Prob. of firing')

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
print(fig1, ['figs' filesep 'Figure1_STA' ], '-dpng', '-r600');
stupid_ratio = 15/16;
myfiguresize = [left, bottom, width*stupid_ratio, height*stupid_ratio];
set(fig1, 'PaperPosition', myfiguresize);
print(fig1, ['figs' filesep 'Figure1_STA'], '-dsvg', '-r600');
