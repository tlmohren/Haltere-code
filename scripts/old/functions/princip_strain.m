function [Ep ] = princip_strain( E,threshold )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% E = [ex,exy,exz,eyy,eyz,ezz];

E = E.*( abs(E)>=threshold);

T = [E(1),  E(2),E(3);...         % create strain tensor
E(2),E(4),E(5);...
E(3),E(5),E(6)];      

[Ve,Ei] = eig(T);                       % determine eigenvalues and eigenvectors 

for kk = 1:3                            % flip eigenvector so largest direction is positive 
%     [~ ,ind]= max(abs(Ve(:,kk)));
%    if Ve(ind,kk)<0
%       Ve(:,kk) = Ve(:,kk).*-1; 
%    end
   
   
   % or use second condition:
   if Ve(1,kk)<=0
      Ve(:,kk) = Ve(:,kk).*-1; 
   end 
   
end

e_vals = diag(Ei);                       % diagonalize
[~ ,ind] = sort(abs(e_vals));           % sort eigenvalues
ind2 = flipud(ind);                     % from large to small
eps = e_vals(ind2);              
eds = reshape(Ve(:,ind2'),1,9);  % flip eigenvecotrs


Ep = [eps(1),eds(1),eds(2),eds(3),eps(2),eds(4),eds(5),eds(6),eps(3),eds(7),eds(8),eds(9)];

end

