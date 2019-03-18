% Matlab function that computes the spacing between the conductors of OHL 
% cables 
%  
%               Programmer: Hemanth Kumar V, Michigan Technological Univ
%               Advisor: Dr Sumit Paudyal, MTU
%               Last Modified: 29th Jan 2015
%

function [Dis,Sij] = Spacingalg(SpacD)

% Importdata.mat file contains the data of configuration, spacing details
% etc
%load importdata.mat

SpacNo = length(SpacD(:,1));

for ii = 1:SpacNo
    %Dis(:,:,SpacD(ii,2)) = zeros(6,6);
    %Sij(:,:,SpacD(ii,2)) = zeros(6,6);
    
    for jj = 1:SpacD(ii,3) % preparing d1, d2, d3 for conductors
        dabcn(1,1,jj) = SpacD(ii,2*(jj)+2)+1i*SpacD(ii,2*(jj)+3);  
    end

    % loop that forms the d matrix and assigns in multidimension array wrt
    % spacing no
    for iii = 1:SpacD(ii,3)
        for jjj = 1:SpacD(ii,3)
            Dis(iii,jjj,SpacD(ii,2)) = abs(dabcn(1,1,iii)-dabcn(1,1,jjj));
        end
    end
    
    % loop that forms the S matrix and assigns in multidimension array wrt
    % spacing no
    for iii = 1:SpacD(ii,3)
        for jjj = 1:SpacD(ii,3)
            Sij(iii,jjj,SpacD(ii,2)) = abs(dabcn(1,1,iii)-conj(dabcn(1,1,jjj)));
        end
    end  
         
end
