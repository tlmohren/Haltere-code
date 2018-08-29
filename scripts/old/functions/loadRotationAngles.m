function [ phi,theta ] = loadExcitationAngles( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    if strcmp(full_dir(end-3:end),'.txt')
        X = importdata( full_dir, 't');
    else

        X = importdata( [full_dir,'.txt'] , 't');
    end

%     first_line = char(X(1));
%     header = all(first_line(3:7) == 'Model');
% 
%     testline = char(X(7));
%     disps = {'Total','X comp','Y comp', ' Z comp'};
%     
%     for j = 1:length(disps)
%         matches = strfind(testline, disps{j});
%         if matches >=0
%             binar(j) = matches;
%         else 
%             binar(j) = 0;
%         end
%     end
% %     binar
%     ndisps = sum(binar>0);
%     c = find(binar);
%     [~,b]=sort(binar(c));
% %     
%     for kk = 1:length(b)
% %         b(kk)
%         order{kk} = disps{c(b(kk))};
%     end
%     
%     if header
%         n_points = length(X)-9;
%         t_steps =  (  length(strsplit(char(X(10))))-3 ) / sum(binar>0);
%         
%         
%         D =zeros(t_steps,n_points,ndisps);
%         
%         for j = 1:n_points
%             y11 = char( X(j+9));
%             y12 = strsplit(y11);
%             z(1,j) = str2double(y12{1});
%             z(2,j) = str2double(y12{2}); 
%             z(3,j) = str2double(y12{3});
% 
%             for k = 1:t_steps
%              
%                 for m = 1:ndisps
% %                     D(k,j,m) = str2double(y12{   3+  k*ndisps + m -ndisps      }); 
%                     D(k,j,m) = str2double(y12{3+    k*ndisps + m -ndisps      }); 
%                 end
%             end
%         end
%     else
%         display('Table not right format?')
%     end
% 
% %     z = round(z,6);
end

