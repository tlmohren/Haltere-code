
% mesh locations
theta = linspace(0,2*pi,17);
theta(1) = 2*pi;
xDes = [0:150:4800];

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
