% Matlab code for calculating the current for capacitirs and loads, it is 
% a part of Three Phase load flow Program 
%               Programmer: Hemanth Kumar V, Michigan Technological Univ
%               Advisor: Dr Sumit Paudyal, MTU
%               Last Modified: 29th Jan 2015



function [Ibus,IL, ILpq, ILz, ILc] = LoadnCap(v,iter,IL, LoadDyn)
load InputFileGenerator.mat
load Inputdata.mat
% load InputFileGenerator.mat
LoadsD = LoadDyn;
% load Inputdata.mat
StudyInfo = input.data.StudyCase;
BusD = input.data.Nodes;
BrchD = input.data.Branch; % Branch Data imports
CapsD = input.data.Caps; % Shunt capacitor details import

DI = [1 0 -1; -1 1 0; 0 -1 1];
D = [1 -1 0; 0 1 -1; -1 0 1];
W = 1/3*[2 1 0; 0 2 1; 1 0 2];

            for ii = 1:length(find(LoadsD(:,1)))  
             ILpq(1:3,1:1,LoadsD(ii,1)) = 0;
             ILz = ILpq;
             ILc = ILpq;
            end


% clear S
 for ii = 1:length(LoadsD(:,1))
     
     if (LoadsD(ii,13)==0)&&(LoadsD(ii,14)==1) % if load data type 
         % col index 13 - '0' - kW & kVAR, '1' - kW& pf, '2' - kVA & pf.
         % col index 14 - '0' - out of status , '1' - in service.
     S(:,:,LoadsD(ii,1)) = [(LoadsD(ii,4) + i*LoadsD(ii,5))*LoadsD(ii,10);
                            (LoadsD(ii,6) + i*LoadsD(ii,7))*LoadsD(ii,11);
                            (LoadsD(ii,8) + i*LoadsD(ii,9))*LoadsD(ii,12)];
     
     elseif (LoadsD(ii,13)==1)&&(LoadsD(ii,14)==1) % kW and pf
     S(:,:,LoadsD(ii,1)) = [(LoadsD(ii,4) + i*(LoadsD(ii,4)/LoadsD(ii,5))*sin(acos(LoadsD(ii,5))))*LoadsD(ii,10);
                            (LoadsD(ii,6) + i*(LoadsD(ii,6)/LoadsD(ii,7))*sin(acos(LoadsD(ii,7))))*LoadsD(ii,11);
                            (LoadsD(ii,8) + i*(LoadsD(ii,8)/LoadsD(ii,9))*sin(acos(LoadsD(ii,9))))*LoadsD(ii,12)];
    
     elseif (LoadsD(ii,13)==2)&&(LoadsD(ii,14)==1) % kVA and pf
     S(:,:,LoadsD(ii,1)) = [ (LoadsD(ii,4)*LoadsD(ii,5)+ i*(LoadsD(ii,4)*sin(acos(LoadsD(ii,5)))))*LoadsD(ii,10);
                             (LoadsD(ii,6)*LoadsD(ii,7)+ i*(LoadsD(ii,6)*sin(acos(LoadsD(ii,7)))))*LoadsD(ii,11);
                             (LoadsD(ii,8)*LoadsD(ii,9)+ i*(LoadsD(ii,8)*sin(acos(LoadsD(ii,9)))))*LoadsD(ii,12)];  
     end
     S(isnan(S))=0;
 end
    
     for ii = 1:length(LoadsD(:,1))  % computes the load current for all loads 
           
         if LoadsD(ii,3)== 0
                f_ziP = 1;
                f_zIp = 0;
                f_Zip = 0;
            elseif LoadsD(ii,3)== 1
                f_ziP = 0;
                f_zIp = 1;
                f_Zip = 0;
            elseif LoadsD(ii,3)== 2
                f_ziP = 0;
                f_zIp = 0;
                f_Zip = 1;
            else
                msgbox(' Yet to model combined load model')
         end
         
            
 %% WYE Connected Load                       
 if (LoadsD(ii,14)==1)&&(LoadsD(ii,2)==1) % load in service and WYE type     
 
     % Calculation of current for PQ type load
for rrr = 1:3
     
    if S(rrr,1,LoadsD(ii,1))~=0
%         ILpq(rrr,1:1,LoadsD(ii,1)) = (conj(S(rrr,:,LoadsD(ii,1))./v(rrr,:,LoadsD(ii,1)))).*(f_ziP*1000);
          ILpq(rrr,1:1,LoadsD(ii,1)) = (conj(S(rrr,:,LoadsD(ii,1))/v(rrr,:,LoadsD(ii,1))))*(f_ziP*1000);

    end
 
    % Calculation of current for constant impedance type load calculation
    % of impedance only during first iteration
if iter == 0   
%     LoadsD(ii,1)
    if S(rrr,:,LoadsD(ii,1))~=0 && StudyInfo(1,5)~=0
%         abs(v(rrr,:,LoadsD(ii,1)))
%         abs(v(rrr,:,LoadsD(ii,1)))*sqrt(3)
        Z(rrr,1,LoadsD(ii,1)) = abs(v(rrr,:,LoadsD(ii,1))).*abs(v(rrr,:,LoadsD(ii,1)))./(conj(S(rrr,:,LoadsD(ii,1)))*1000);

%         ILM(rrr,1,LoadsD(ii,1)) = abs(conj(S(rrr,:,LoadsD(ii,1))...
%                         *1000./v(rrr,:,LoadsD(ii,1)))).*f_zIp;
                    
        ILM(rrr,1,LoadsD(ii,1)) = conj(S(rrr,:,LoadsD(ii,1))*1000/v(rrr,:,LoadsD(ii,1)))*f_zIp;
    end            
             
    if S(rrr,:,LoadsD(ii,1))~=0 && StudyInfo(1,5)==0
        [r0,c0,v0] = find(BusD(:,1) == LoadsD(ii,1));
%          LoadsD(ii,1)
        vi = BusD(r0,4)/sqrt(3);
        Z(rrr,1,LoadsD(ii,1)) = abs(vi*vi*10^6)/(conj(S(rrr,:,LoadsD(ii,1))*1000));
        ILM(rrr,1,LoadsD(ii,1)) = abs(conj(S(rrr,:,LoadsD(ii,1))...
                                            *1000./(vi*1000))).*f_zIp;
%         ILM(rrr,1,LoadsD(ii,1)) = conj(S(rrr,:,LoadsD(ii,1))*1000/(vi*1000))*f_zIp;
    end 
end
     
    



    if S(rrr,:,LoadsD(ii,1))~=0
        ILz(rrr,1:1,LoadsD(ii,1)) = (v(rrr,:,LoadsD(ii,1))/Z(rrr,:,LoadsD(ii,1))).*f_Zip;
    end
  
% Calculation of current for constant current type load  
if S(rrr,:,LoadsD(ii,1))~=0
    ILc(rrr,1:1,LoadsD(ii,1)) =  ILM(rrr,1,LoadsD(ii,1))...
    *cosd((angle(v(rrr,:,LoadsD(ii,1)))*180/pi)-(angle(S(rrr,:,LoadsD(ii,1)))*180/pi))+i*ILM(rrr,1,LoadsD(ii,1))*sind((angle(v(rrr,:,LoadsD(ii,1)))*180/pi)-(angle(S(rrr,:,LoadsD(ii,1)))*180/pi));
end
                           
 end
% Total Load current


IL(:,:,LoadsD(ii,1)) = ILpq(1:3,1:1,LoadsD(ii,1))+ ILz(1:3,1:1,LoadsD(ii,1))+ILc(1:3,1:1,LoadsD(ii,1));
                        
 end
 
         
 %% Delta connected Load       
if (LoadsD(ii,2)== 0)&&(LoadsD(ii,14)==1)  % if the load == Delta connection & in status

 % Calculation of current for PQ type load 
 Vde(1:3,1:1,LoadsD(ii,1)) = D*v(1:3,1:1,LoadsD(ii,1));
 
 for rrr=1:3
     
     if S(rrr,:,LoadsD(ii,1))~=0 
        ILpq(rrr,1,LoadsD(ii,1)) = (conj(S(rrr,1,LoadsD(ii,1))./(Vde(rrr,1,LoadsD(ii,1))))).*(f_ziP*1000);
     end
 
    % Calculation of current for constant impedance type load         
    % calculation of impedance only during first iteration  
    if iter == 0        
       
         if S(rrr,:,LoadsD(ii,1))~=0 
    Z(rrr,1:1,LoadsD(ii,1)) = abs(Vde(rrr,:,LoadsD(ii,1)).*Vde(rrr,:,LoadsD(ii,1)))...
                                    ./(conj(S(rrr,:,LoadsD(ii,1)))*1000);
    ILM(rrr,1:1,LoadsD(ii,1)) = abs(conj(S(rrr,:,LoadsD(ii,1))...
                                    *1000./(Vde(rrr,1,LoadsD(ii,1))))).*f_zIp;
         end
             
             
         if S(rrr,:,LoadsD(ii,1))~=0 && StudyInfo(1,5)==0
            [r0,c0,v0] = find(BusD(:,1) == LoadsD(ii,1));
            vi = BusD(r0,4);
            Z(rrr,1,LoadsD(ii,1)) = abs(vi*vi*10^6)./(conj(S(rrr,:,LoadsD(ii,1)))*1000);
          ILM(rrr,1,LoadsD(ii,1)) = abs(conj(S(rrr,:,LoadsD(ii,1))*1000./(vi*1000))).*f_zIp;

%             ILM(rrr,1,LoadsD(ii,1)) = conj(S(rrr,:,LoadsD(ii,1))*1000/(vi*1000))*f_zIp;
         end            
    end
 
  if S(rrr,:,LoadsD(ii,1))~=0
    ILz(rrr,1:1,LoadsD(ii,1)) = ((Vde(rrr,1,LoadsD(ii,1)))/Z(rrr,1,LoadsD(ii,1)))*f_Zip;
  end
  
% Calculation of current for constant current type load  
     if S(rrr,:,LoadsD(ii,1))~=0
        ILc(rrr,1,LoadsD(ii,1)) = ILM(rrr,1,LoadsD(ii,1))*cosd((angle(Vde(rrr,:,LoadsD(ii,1)))*180/pi)-(angle(S(rrr,:,LoadsD(ii,1)))*180/pi))+i*ILM(rrr,1,LoadsD(ii,1))*sind((angle(Vde(rrr,:,LoadsD(ii,1)))*180/pi)-(angle(S(rrr,:,LoadsD(ii,1)))*180/pi));
     end
 end                        

% Total Load current
IL(:,:,LoadsD(ii,1)) = DI*ILpq(1:3,1:1,LoadsD(ii,1))+ DI*ILz(1:3,1:1,LoadsD(ii,1))+DI*ILc(1:3,1:1,LoadsD(ii,1)) ;
                                             
end % end of delta load current
  
end                                       
iter = iter+1;                                        

%% Shunt Capacitor modelling
% for cc = 1:length(find(CapsD(:,1)))     
% % Shunt Capacitor connected in WYE Mode
% % if CapsD(cc,2)==1 && CapsD(cc,11)~=0 % wye connected
% %     % Bc = kVar/(kV_LN^2 * 1000)
% %       Bc(1,1) = CapsD(cc,6)/(power(CapsD(cc,4),2)*1000) 
% %       Bc(2,1) = CapsD(cc,8)/(power(CapsD(cc,4),2)*1000)
% %       Bc(3,1) = CapsD(cc,10)/(power(CapsD(cc,4),2)*1000)
% %       
% %     Vc = v(:,:,CapsD(cc,1))  
% %     IC(:,:,CapsD(cc,1)) = [i*Bc(1,1)*Vc(1,1);i*Bc(2,1)*Vc(2,1);i*Bc(3,1)*Vc(3,1)]   ;
% % end
% 
% % % Shunt Capacitor connected in Delta Mode
% % if CapsD(cc,2) == 0 && CapsD(cc,11)~=0  % 0-delta connected
% %   %Warning ! the kV rating of the delta connected cap bank should be line to line voltage only
% %     % Bc = kVar/(kV_LL^2 * 1000)
% %       Bc(1,1) = CapsD(cc,6)/(power(CapsD(cc,4),2)*1000); 
% %       Bc(2,1) = CapsD(cc,8)/(power(CapsD(cc,4),2)*1000);
% %       Bc(3,1) = CapsD(cc,10)/(power(CapsD(cc,4),2)*1000);
% %  
% % Vc = v(:,:,CapsD(cc,1)); % Check if the voltage obtained is the 'Line to Line' or the 'Line to Neutral' voltgaes  
% % IC(:,:,CapsD(cc,1)) = DI*[i*Bc(1,1)*Vc(1,1);i*Bc(2,1)*Vc(2,1);i*Bc(3,1)*Vc(3,1)];
% %   
% % end 
%  
% end

%% Summing all the shunt loads and capacitor current into a single constant called I_bus
for ii = 1:length(BusD(:,1))
Ibus(:,:,BusD(ii,1)) = IL(:,:,BusD(ii,1));%+IC(:,:,BusD(ii,1));
% The load current and the capacitor currents are added together at the node
end
end
