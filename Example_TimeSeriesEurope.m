clc
clear all
DI = [1 0 -1; -1 1 0; 0 -1 1];
D = [1 -1 0; 0 1 -1; -1 0 1];
W = 1/3*[2 1 0; 0 2 1; 1 0 2];

% This particular piece of code executes the importing data and information 
% run InputFileGeneratorR8.m   



% SCase information indicates which Case Scenario to opt for in Studyinfo 
% data. For current version it is always 1 
SCase = 1; 
% time_end = 2; % 1 indicates the first snapshot i.e., data available in LoadsD is considered.
              % not equal 1 indicates the load time-series data is
              % avialable. But the load data needs to be in .mat file
              % format. This will be loaded in threephase function
              
plot_fig = 1; % Make this 1 if plot function below needs to be executed !


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
 hr=1:5 ; 
time_end = length(hr);

[TreeTab] = TreeAlgR6(SCase);


% Calling on to three-phase load flow function
    disp(hr);
    [v,I, IL, ILpq, ILc, ILz] = ThreePhLF(hr, SCase, CapsD, RegD);
  %%  
  
  for hh=1:time_end
    for ii = 1:length(BusD(:,1))
        vll(:,:,BusD(ii,1),hr(hh)) = D*v(:,:,BusD(ii,1));
        vln(:,:, BusD(ii,1),hr(hh)) = v(:,:,BusD(ii,1));
    end
    
    for ii = 1:length(BrchD(:,1))
        Ibr(:,:, BrchD(ii,1), hr(hh)) = I(:,:,BrchD(ii,1),1);
    end         
  end


% % % Calling on to three-phase load flow function
% % for hr = 1:time_end
% %     disp(hr)
% %     [v,I, IL, ILpq, ILc, ILz] = ThreePhLF(hr, SCase, CapsD, RegD);
% %     for ii = 1:length(BusD(:,1))
% %         vll(:,:,BusD(ii,1),hr) = D*v(:,:,BusD(ii,1));
% %         vln(:,:, BusD(ii,1),hr) = v(:,:,BusD(ii,1));
% %     end
% %     
% %     for ii = 1:length(BrchD(:,1))
% %         Ibr(:,:, BrchD(ii,1), hr) = I(:,:,BrchD(ii,1),1);
% %     end          
% % end




% [v,I, IL, ILpq, ILc, ILz] = ThreePhLF(hr, SCase, CapsD, RegD);
 
% run OutputFileGenerator.m

%

%% Ploting
% the hr variable here needs to be entered based on the user input above.
if plot_fig == 1
    
        Ia(1, hr) = abs(Ibr(1,1, 1, hr));
        Ib(1, hr) = abs(Ibr(2,1, 1, hr));
        Ic(1, hr) = abs(Ibr(3,1, 1, hr));

        %  It is upto user to plot the no of attribues based on need. I havenot generalized and kept it open
        plot(hr, Ia(1, hr),'r', hr, Ib(1, hr),'b', hr, Ic(1, hr), 'g');
        xlabel('time');
        ylabel('Amperes');
        % set(gca,'Color','k')
        grid on
        
        
        figure
        V1(1, hr) = abs(vln(1,1, LoadsD(1,1), hr));
        V32(1, hr) = abs(vln(3,1, LoadsD(32,1), hr));
        V53(1, hr) = abs(vln(2,1, LoadsD(53,1), hr));
        V1ang(1, hr) = angle(vln(1,1, LoadsD(1,1), hr))*180/pi;
        V32ang(1, hr) = angle(vln(3,1, LoadsD(32,1), hr))*180/pi;
        V53ang(1, hr) = angle(vln(2,1, LoadsD(53,1), hr))*180/pi;
        
        %  It is upto user to plot the no of attribues based on need. I havenot generalized and kept it open
        plot(hr, V1(1, hr),'r', hr, V53(1, hr),'b', hr, V32(1, hr), 'g');
        ylim([235 255])
        legend('Load 1(ph A)','Load 53(ph B)', 'Load 32(ph C)')
        ylabel('Load voltage')
        xlabel('Time (mins)')
end
%%
save ('timeseries_eurp.mat')
% % %% Computing Substation Real and Reactive power through Branch-1
% % hr = 1:time_end
%% 
for hh=1:time_end
    vnull(:,:) = W*vll(:,:,1,hh);
    inull(:,:) = conj(Ibr(:,:,1,hh));
    Sgrid(1,1, hh) = vnull(1,1)*inull(1,1);
    Sgrid(2,1, hh) = vnull(2,1)*inull(2,1);
    Sgrid(3,1, hh) = vnull(3,1)*inull(3,1);
    Pa(1,hh) = real(Sgrid(1,1, hh));
     Pb(1,hh) = real(Sgrid(2,1, hh));
      Pc(1,hh) = real(Sgrid(3,1, hh));
      Ptotal(1,hh) = Pa(1,hh)+Pb(1,hh)+Pc(1,hh);
      Qtotal(1,hh) = imag(Sgrid(1,1, hh))+imag(Sgrid(2,1, hh))+imag(Sgrid(3,1, hh));
end
 
%%  Substation Active and reactive power 
if plot_fig == 1
        plot(hr, Pa,'r', hr, Pb,'b', hr, Pc, 'g')
        xlabel('time')
        ylabel('kW')
        % set(gca,'Color','k')
        grid on

        figure
        plot(hr, Ptotal,'r', hr, Qtotal, 'b')
        xlabel('Time')
        ylabel('kW/kVAr')
        grid on

        legend('Real Power', 'Reactive Power','Location','northwest')
end