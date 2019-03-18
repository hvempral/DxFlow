clc
clear 

% sequence to Zabc
Zpos= 0.322+0.074i






Zzero =	0.804+0.093i
Zs = (Zpos*2 +Zzero)/3;
Zm = Zs-Zpos;

% Zabc = [Zs Zm Zm;Zm Zs Zm; Zm Zm Zs]

% Zabc1 = [0.4576+1.078i 0.1559+0.5017i 0.1535+0.3849i;
%          0.1559+0.5017i 0.4666+1.0482i 0.158+0.4236i;
%          0.1535+0.3849i 0.158+0.4236i 0.4615+1.0651i]
Z012 = [Zzero 0 0; 0 Zpos 0; 0 0 Zpos]     
     
a = -0.5 + 0.866i;
A = 1/3*[1 1 1; 
         1 a a*a; 
         1 a*a a]
     
%      Z012 = A*Zabc*inv(A)
%      Zabc
     
%  Z = [real(Zabc(1,1)), imag(Zabc(1,1)) ,real(Zabc(1,2)), imag(Zabc(1,2)) ,real(Zabc(1,3)), imag(Zabc(1,3)),real(Zabc(2,2)), imag(Zabc(2,2)),real(Zabc(2,3)), imag(Zabc(2,3)) ,real(Zabc(3,3)), imag(Zabc(3,3))]
     
     
 
 As = [1 1 1; 
         1 a*a a; 
         1 a a*a]
     
     Zabc = As*Z012*inv(As)
     
     Z = [real(Zabc(1,1)), imag(Zabc(1,1)) ,real(Zabc(1,2)), imag(Zabc(1,2)) ,real(Zabc(1,3)), imag(Zabc(1,3)),real(Zabc(2,2)), imag(Zabc(2,2)),real(Zabc(2,3)), imag(Zabc(2,3)) ,real(Zabc(3,3)), imag(Zabc(3,3))]
