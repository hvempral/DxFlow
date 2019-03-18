function [Ibus, IC] = Caps(v,IC,IL,Ibus, CapsD, BusD)
%% Shunt Capacitor modelling

DI = [1 0 -1; -1 1 0; 0 -1 1];
D = [1 -1 0; 0 1 -1; -1 0 1];
W = 1/3*[2 1 0; 0 2 1; 1 0 2];

for cc = 1:length(find(CapsD(:,1)))     
%% Shunt Capacitor connected in WYE Mode
if CapsD(cc,11)~=0
    if CapsD(cc,12)==0 % 0-Fixed capacitor, rated kVAR in CapsD(6,8,10) is used
            % Bc = kVar/(kV_LN^2 * 1000)
    Bc(1,1) = CapsD(cc,6)/(power(CapsD(cc,4),2)*1000);
    Bc(2,1) = CapsD(cc,8)/(power(CapsD(cc,4),2)*1000);
    Bc(3,1) = CapsD(cc,10)/(power(CapsD(cc,4),2)*1000);
    else
    Bc(1,1) =  CapsD(cc,13)*CapsD(cc,14)/(power(CapsD(cc,4),2)*1000);
    Bc(2,1) =  CapsD(cc,13)*CapsD(cc,15)/(power(CapsD(cc,4),2)*1000);
    Bc(3,1) =  CapsD(cc,13)*CapsD(cc,16)/(power(CapsD(cc,4),2)*1000);        
    end
    
    Vc = v(:,:,CapsD(cc,1)); 
    % Check if the voltage obtained is the 'Line to Line' or the 
    % 'Line to Neutral' voltgaes  

    if CapsD(cc,2)==1  % wye connected
        IC(:,:,CapsD(cc,1)) = [ i*Bc(1,1)*Vc(1,1);
                                i*Bc(2,1)*Vc(2,1);
                                i*Bc(3,1)*Vc(3,1)] ;                 
    else   % Shunt Capacitor connected in Delta Mode
    %Warning ! The kV rating of the delta connected cap bank 
    %          should be line to line voltage only
    % Bc = kVar/(kV_LL^2 * 1000)
    IC(:,:,CapsD(cc,1)) = DI*[  i*Bc(1,1)*Vc(1,1);
                                i*Bc(2,1)*Vc(2,1);
                                i*Bc(3,1)*Vc(3,1)]; 
    end 
 
end

end

    for ii = 1:length(BusD(:,1))
    Ibus(:,:,BusD(ii,1)) = IL(:,:,BusD(ii,1))+IC(:,:,BusD(ii,1));
    % The load current and the capacitor currents are added together at the node
    end
end % end of function
