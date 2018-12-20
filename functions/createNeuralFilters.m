function [ STAFunc,NLDFunc] = createNeuralFilters( STAfreq,STAwidth,STAdelay,NLDgrad,NLDshift )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    STAFunc = @(t) -cos( STAfreq*(t+STAdelay)  ).*exp(-(t+STAdelay).^2 / STAwidth^2);
    NLDFunc = @(s) ( 1./ (1+ exp(-NLDgrad.*(s-NLDshift)) ) )  ; 

%     STAt = -39:0.1:0;   
%     NLDrange = -1:0.01:1;
%     figure(100);
%     subplot(211)
%         plot(STAt, STAFunc(STAt) );hold on;drawnow
%     subplot(212)
%         plot(NLDrange, NLDFunc(NLDrange) );hold on;drawnow; grid on
end

