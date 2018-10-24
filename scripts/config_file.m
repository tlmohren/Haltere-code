
% mesh locations
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;
xDes = [0:150:4800];


fsz = 7;
set(0,'DefaultAxesFontSize',fsz)% .
set(0,'DefaultLegendFontSize',fsz)% .
% set(0,   'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
%       'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here

% Figure text font multiplier 
set(0,'DefaultAxesLabelFontSize', 1)
set(0,'DefaultAxesTitleFontSize', 1)


circleDistance = 3000;               % distance from base to haltere 
circleRadius = 150;                 % radius of haltere   
        
n = 16;
% colorscheme
redPurple = [158,1,66
213,62,79
244,109,67
253,174,97
254,224,139
255,255,191
230,245,152
171,221,164
102,194,165
50,136,189
94,79,16];
% 
strainScheme = colorSchemeInterp(redPurple/255, 500);

% coloring of plot
surfParamBackground = {'FaceAlpha',1,'EdgeAlpha',0.4};
surfParamBackgroundTwistStalk = {'FaceAlpha',1,'EdgeAlpha',0.4};
surfParamForeground = {'FaceAlpha',1,'EdgeAlpha',0.4};

% axis opts 

% used indices for plotting 
% It = 1:160;
len = 101;
start = 35;
It = start:(start+len-1);
t_plot = (0:len-1)*0.001;

% axOpts = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15]}; 
% axOptsNoT = {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15],'XTickLabel',{'','','',''} }; 

% t = 0:0.001:0.35;

deformLabels = {'$\Delta \phi$','$\Delta \theta$','$\Delta \gamma$'};
        

% axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.10],'XTickLabel',{'','','',''} ,...
%                'YLim',[-1,1]*0.005}; 
% axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.10],'XTickLabel',{'','','',''} ,...
%                'YLim',[-1,1]*2e-4}; 
% axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.10] ,...
%                'YLim',[-1,1]*2e-6};           
           
           
axOpts_dphi= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15],'XTickLabel',{'','','',''} ,...
               'YLim',[-1,1]*0.005}; 
axOpts_dtheta= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15],'XTickLabel',{'','','',''} ,...
               'YLim',[-1,1]*3e-4}; 
axOpts_dgamma= {'XGrid','On','XLim',[0,t_plot(end)],'XTick',[0:0.05:0.15] ,...
               'YLim',[-1,1]*3e-6}; 
           
           
           