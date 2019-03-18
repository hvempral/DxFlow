clc
clear all

DI = [1 0 -1; -1 1 0; 0 -1 1];
D = [1 -1 0; 0 1 -1; -1 0 1];
W = 1/3*[2 1 0; 0 2 1; 1 0 2];

% This particular piece of code executes the importing data and information 
run InputFileGeneratorR8.m   

% SCase information indicates which Case Scenario to opt for in Studyinfo 
% data. For current version it is always 1 
SCase = 1; 
[TreeTab] = TreeAlgR6(SCase);
% Loading basic element info.              
load Inputdata.mat 
BusD = input.data.Nodes;
BrchD = input.data.Branch; % Branch Data imports
ConfgD = input.data.Configuration ;% Branch data configuration import
LoadsD = input.data.Loads; % Load data import
TrfD = input.data.Transformers; % Transformer data import
StudyInfo = input.data.StudyCase; % STudy cases imports
% Can participates in optimization
CapsD = input.data.Caps; % Capacitor data 
RegD = input.data.Regulator; % Regulator data 

% this is used for time-series. Can be initialized as 1 for optimization case
hr=1; 
% Calling on to three-phase load flow function
tic
[v,I, IL, ILpq, ILc, ILz] = ThreePhLF(hr, SCase, CapsD, RegD);
toc    


for ii =1:size(BusD(:,1))
Aph(ii,1) = abs(v(1,1,BusD(ii,1)))/(BusD(ii,4)*1000);
end

plot(Aph,'MarkerSize',3,'Marker','o','Color','r', 'markerfacecolor', 'b')
xlim([0 740])
ylim([0.9 1.05])
xlabel('Node number')
ylabel('Node voltages in pu')
grid on

% % Total load calculated at each load bus
% % real and reactive power at each of the load type and load model is
% % accounted here (tested)
% totalload = zeros(3,1);
% for ii=1:length(LoadsD(:,1))
%     Iload = ILpq(:,1, LoadsD(ii,1))+ILz(:,1, LoadsD(ii,1))+ILc(:,1, LoadsD(ii,1));
%     vloadll = D*v(:,:,LoadsD(ii,1));    
%     if LoadsD(ii,2)==0
%         totalload = totalload + vloadll.*conj(Iload);
%     else 
%         totalload = totalload + v(:,:, LoadsD(ii,1)).*conj(Iload);
%     end
%     
%     totalload;
% end
% 
% %% Calculation of branch feeding system from source 
% srcbus = 60
% Ssource = v(:,:, srcbus).*conj(Ibr(:,:, 94,hr))
% TotalLoss = (real(Ssource - totalload)) % sum of each phase real powers (W)
% % run OutputFileGenerator.m

