% Matlab code for voltage regulator function created 10/3/2014, as a part of 
% Three Phase load flow Program 
%               Programmer: Hemanth Kumar V, Michigan Technological Univ
%               Advisor: Dr Sumit Paudyal, MTU
%               Last Modified: 29th Jan 2015


 function [Tap, Brchs] = VoltReg(RegD, BusD, Brchs, iter) % ,v, Ibus, T, I)

Vrel = [];
RegNo = length(find(RegD(:,1)));
RegId = RegD(:,1);

if RegNo ~= 0
    BW = RegD(:,9);
    Pt = RegD(:,10);
    Ct = RegD(:,11);
    Vmin = RegD(:,12);
    Vmax = RegD(:,13);
    NoStp = RegD(:,14);

for reii = 1:RegNo

    PerStp(reii,1) = ((Vmax(reii,1)-Vmin(reii,1))/NoStp(reii,1))*100;
    [r,c1,xxx] = find(BusD(:,1)==RegD(reii,4));
    StepV = BusD(r,4)*1000/(sqrt(3)*Pt(reii,1));
    PerStepV = PerStp(reii,1)*StepV/100;
    
    %% for first iter and fixed tap settings
if iter==0 % For first iter it is reqd to 
%         Compute the abcd parameters with the given taps and also if 
%         the settings are forced this part of code has to run
%         For RegD(:,18)
%         0- Tap settings are forced
%         1- Tap settings are free to change in  load flow    2- tap settings are forced but allow during
%         optimization 

    Tap(:,:,reii) = [ RegD(reii, 15); 
                      RegD(reii, 16); 
                      RegD(reii, 17)];
  
    if RegD(reii, 6)== 1 % 0- WYE type A & 1- Delta Type B
%         r_a(:,:,RegId(reii,1)) = [1-PerStp(reii,1)/100*Tap(1,1,reii)     0                               0;
%                          0                          1-PerStp(reii,1)/100*Tap(2,1,reii)          0; 
%                          0                                      0                       1-PerStp(reii,1)/100*Tap(3,1,reii)]
        r_a(:,:,RegId(reii,1)) = [1-0.005888*Tap(1,1,reii)     0                               0;
                         0                          1-0.005888*Tap(2,1,reii)          0; 
                         0                                      0                       1-0.005888*Tap(3,1,reii)];
 
                     
                     r_b(:,:,RegId(reii,1)) = zeros(3,3);
        r_c(:,:,RegId(reii,1)) = zeros(3,3);
        
        
        r_d(:,:,RegId(reii,1)) = [ 1/r_a(1,1,RegId(reii,1))             0           0;
                           0                   1/r_a(2,2,RegId(reii,1))    0;
                           0                   0                1/r_a(3,3,RegId(reii,1))];
        
        
        r_d(~isfinite(r_d))=0;
        
        r_A(:,:,RegId(reii,1)) =  r_d(:,:,RegId(reii,1));
        r_B(:,:,RegId(reii,1)) = zeros(3,3);
%         test(reii,1) = 10              
    else
        disp(' yet to model the closed delta confg (Type B)')
    end
   
  
    elseif RegD(reii,18)==0 % If Tap settings are forced
        Tap(:,:,reii) = [RegD(reii, 15); RegD(reii, 16); RegD(reii, 17)]
        
        if RegD(reii, 6)== 1 % 0- WYE type A & 1- Delta Type B
            r_a(:,:,RegId(reii,1)) = [1-0.005888*Tap(1,1,reii) 0 0; 
                0 1-0.005888*Tap(2,1,reii) 0; 
                0 0 1-0.005888*Tap(3,1,reii)];
            r_b(:,:,RegId(reii,1)) = zeros(3,3);
            r_c(:,:,RegId(reii,1)) = zeros(3,3);
            r_d(:,:,RegId(reii,1)) = inv(r_a(:,:,RegId(reii,1)));
            r_A(:,:,RegId(reii,1)) = inv(r_a(:,:,RegId(reii,1)));
            r_B(:,:,RegId(reii,1)) = zeros(3,3);
%             test(reii,1) = 10
        else
            disp(' yet to model the closed delta confg (Type B)')
        end
    
        %% for changing the tap during load flow and short circuit given that 
    % iter is not equal to zero
elseif ((iter~=0)&&(RegD(reii,18)==1))%//(RegD(rii,18)==2))
        % Tap settings are free to change, which is altered only during the
        % load flow
    Vpt = v(:,:,r)./Pt(reii,1)   ;
    Ireg = I(:,:,r,3)/Ct;
    abs(Vpt);
    angle(Vpt);
%     test(reii,1) = 20;
    
    VrgSmin = RegD(reii,28)-RegD(reii,9)/2; % Voltage 
    % regulator setting at lower end
    
    VrgSmax = RegD(reii,28)+RegD(reii,9)/2; % voltage 
    %regulator setting for upper end
    Zcomp = zeros(3,3)
    
for rij = 1:3
    if RegD(reii,18+rij)==0 % rx settings are taken from database
       Rreg = RegD(reii,20+(2*rij))
       Xreg = RegD(reii,20+(2*rij+1))
       Zc_reg = (Rreg+i*Xreg)/RegD(reii,30)
       Zcomp(rij,rij,reii) = Zc_reg           
               
    elseif RegD(reii,18+rij)==1 % Take r x from the line data
       % Rreg = Line Data
       % Xreg = Line Data
       % Zc_reg = (Rreg+i*Xreg)/RegD(Reii,30)
       % Zcomp(:,:,reii) = Zc_reg(rij,rij)   
        
    end
 % Tap setting calculations
    Vrel(rij,1) = Vpt(rij,1)-Zcomp(rij,rij)*Ireg(rij,1)  ;
    abs(VrgSmax);
    abs(Vpt);
    
        if (RegD(reii,28) <= abs(Vrel(rij,1)))
        Tap(rij,1,reii) = -1*round((VrgSmax-abs(Vrel(rij,1)))/PerStepV); %VrgSmax
        else
            RegD(reii,28);
            VrgSmin;
            abs(Vrel(rij,1));
            PerStepV;
            Tap(rij,1,reii) = round((VrgSmin-abs(Vrel(rij,1)))/PerStepV);
       end
end
% Zcomp
% abs(v(:,:,r))
% abs(Vpt)
% abs(Vrel)

    if RegD(reii, 6)==1 % 0- WYE type A & 1- Delta Type B
    r_a(:,:,RegId(reii,1)) = [1-PerStp(reii,1)/100*Tap(1,1,reii) 0 0; 0 1-PerStp(reii,1)/100*Tap(2,1,reii) 0; 0 0 1-PerStp(reii,1)/...
        100*Tap(3,1,reii)];
    r_d(:,:,RegId(reii,1)) = inv(r_a(:,:,RegId(reii,1)));
    r_b(:,:,RegId(reii,1)) = zeros(3,3);
    r_c(:,:,RegId(reii,1)) = zeros(3,3);
    r_A(:,:,RegId(reii,1)) = inv(r_a(:,:,RegId(reii,1)));
    r_B(:,:,RegId(reii,1)) = zeros(3,3);
    else
        disp(' yet to model the closed delta confg (Type B)')
    end  
    
    
    
 %% forcing the tap settings from the optimization code   
elseif RegD(reii,18)==2 % for Tap settings are forced during optimization
%     test(reii,1) = 30;
    if RegD(reii, 6)==0 % 0- WYE type A & 1- Delta Type B
        r_a(:,:,RegId(reii,1)) = [1-PerStp(reii,1)/100*Tap(1,1,reii) 0 0; 0 1-PerStp(reii,1)/100*Tap(2,1,reii) 0; 0 0 1-PerStp(reii,1)/...
            100*Tap(3,1,reii)]
        r_d(:,:,RegId(reii,1)) = inv(r_a(:,:,RegId(reii,1)))
        r_A(:,:,RegId(reii,1)) = inv(r_a(:,:,RegId(reii,1)))
        r_b(:,:,RegId(reii,1)) = zeros(3,3);
        r_c(:,:,RegId(reii,1)) = zeros(3,3);
        r_B(:,:,RegId(reii,1)) = zeros(3,3)
    else
        disp(' yet to model the closed delta confg (Type B)')
    end      
    
    
end
 

end

else
 r_A = []
 r_B = []
 r_a = []
 r_b = []
 r_c = []
 r_d = []
 
 
end

 Brchs.Reg.a = r_a;
  Brchs.Reg.b = r_b;
   Brchs.Reg.c = r_c;
    Brchs.Reg.d = r_d;
     Brchs.Reg.A = r_A;
      Brchs.Reg.B = r_B;
 