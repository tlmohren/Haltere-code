function [Ep,Ed ] = calcPrincipleStrainSeries( E,threshold,n_strains )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% E = [ex,exy,exz,eyy,eyz,ezz];

%     E = E.*( abs(E)>=threshold);
%%% E must have size(time, locations )
    for j = 1:size(E,1)
        Etemp = squeeze(E(j,:));
        T = [Etemp(1),  Etemp(2),Etemp(3);...         % create strain tensor
        Etemp(2),Etemp(4),Etemp(5);...
        Etemp(3),Etemp(5),Etemp(6)];      

        [Ve,Ei] = eig(T);                       % determine eigenvalues and eigenvectors 

        for kk = 1:3                            % flip eigenvector so largest direction is positive 
%             [~ ,ind]= max(abs(Ve(:,kk)));
%            if Ve(ind,kk)<0
%               Ve(:,kk) = Ve(:,kk).*-1; 
%            end
% 
% 
%            % or use second condition:
%            if Ve(1,kk)<=0
%               Ve(:,kk) = Ve(:,kk).*-1; 
%            end 

        end

        e_vals = diag(Ei);                       % diagonalize
        [~ ,ind] = sort(abs(e_vals));           % sort eigenvalues
        ind2 = flipud(ind);                     % from large to small
        eps(j,:) = e_vals(ind2);              
        eds(j,:) = reshape(Ve(:,ind2'),1,9);  % flip eigenvecotrs

    end
%     figure()
%     plot(
%     
    
    %% Correct for switching behavior
%     n_switches = 6;
       der = diff(eps(:,1));
       der2 = diff(eps(:,2));
%         [V,I] = sort(abs(der),'descend');
%         I_sort = sort(I(1:n_switches),'ascend');
%        for j = 1:n_switches
%           eps(I_sort(j)+1:end) = - eps(I_sort(j)+1:end);
%           
%                 figure();
%                    der = diff(eps(:,1));
%                 subplot(211)
%                     plot(eps(:,1))
%                     hold on
%                     
%                     plot(eps(:,2))
%                 subplot(212); hold on 
%                     plot(der)
%                     plot(der2)
%                     der(3:end) 
              
%     for j = 1:size(E,1)
%        der = diff(eps(:,1));


% 




    if n_strains == 1
        Ep = [eps(:,1)];
        Ed= [eds(:,1),eds(:,2),eds(:,3)];
    elseif n_strains == 2
        Ep = [eps(:,2)];
        Ed= [eds(:,4),eds(:,5),eds(:,6)];
%         Ep = [eps(2),eds(4),eds(5),eds(6)];
    else
        Ep = [eps(1),eds(1),eds(2),eds(3),eps(2),eds(4),eds(5),eds(6),eps(3),eds(7),eds(8),eds(9)];
    end
end

