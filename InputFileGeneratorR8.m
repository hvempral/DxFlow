%% Input file generator - Transmission Lines, Transformer

clc
clear
ver = 'new'

%% revision of three phase load flow program
[filename pathname] = uigetfile({'*.xlsx'},'File Selector');
fullpathname = strcat(pathname, filename);
out = strcat(pathname, 'OUT'); % to develop a folder path & link
recycle = strcat(pathname,'BIN')

mkdir(recycle);
mkdir(out);

% Saving path, links, dirrectory names to command.mat file.
save ('command.mat','filename','fullpathname','pathname','out','recycle')

%fullpathname = strcat(pathname, filename);


disp('Importing Section');
% The below line loads DATA from excel sheet to input structure array with
% each sheet in fields
input = importdata(filename)
disp('End of Importing Section');
tser = tic;

%% User input for process
SCase = 1;

% Extracting the data of power system elements
BusD = input.data.Nodes; % Bus Data import
BrchD = input.data.Branch; % Branch Data imports
ConfgD = input.data.Configuration ;% Branch data configuration import
UdbD = input.data.UDB_Line ;
RegD = input.data.Regulator; % Regulator data import
LoadsD = input.data.Loads; % Load data import
TrfD = input.data.Transformers; % Transformer data import
CondLibD = input.data.Conductortype; % conductor type details import
CapsD = input.data.Caps; % Shunt capacitor details import
SpacD = input.data.Spacing; % COnductor spacing details import
StudyInfo = input.data.StudyCase; % STudy cases imports

telp = toc(tser);

if (ver=='new')
UdbT = input.data.UDB_trafo; % Transformer user defined impedance block
end



disp('End of Importing Section');
%save ('importdata.mat','input')

%% Re-modifying the inputfile generator from legacy code of 2015s
freq = StudyInfo(SCase,6);
row =  StudyInfo(SCase,7);
%   end

iter = 0;
%%
W = 1/3*[2 1 0; 
         0 2 1; 
         1 0 2];
     
D = [1 -1  0; 
     0  1 -1; 
    -1  0  1];       

%% Processing the data for finding the series impedance of the Line
stp1 = waitbar(0,'Processing Data......');
waitbar(0.1);
disp(' #################################################################')
disp('                        Series impedance Section                  ');

ConfgId = ConfgD(:,2); % stores the configuration Ids 
% CondLibNo = length(CondLibD(:,1)); 

% BrchNo--> No of Branches , BrchId --> Branch Ids

BrchId = BrchD(:,1);

% Distance algorithm between the conductors and D and S matrix calculation

[Dis,Sij] = Spacingalg(SpacD);
z=zeros(6,6); % pre allocating the matrices for faster performance 

z_abc = zeros(3,3);

for iii = 1:length(ConfgD(:,2))   
    n = ConfgD(iii,11); % the primitive size of the matrix
    Dod = ConfgD(iii,12);% Outside diameter (valid for cable only!!)
    k = ConfgD(iii,13); % No of strands for cable type system 
    
     for ii=1:n % 6 indicates the formation of primitive matrix which 
         % will be different for different scenarion
         % Change with dynamic input for this
           for jj=1:n
                if  ConfgD(iii,10) == 0 % OH Transmission Line        
                     [r2,c2,v2] = find(CondLibD(:,1)== ConfgD(iii,jj+5));             
                     rii =  CondLibD(r2,4);
                     GMRii = CondLibD(r2,6);

                        if (ConfgD(iii,jj+14)==1)&&(ConfgD(iii,ii+14)==1)
                            if (ii==jj)
                                 z(ii,jj,ConfgD(iii,1)) = rii+(0.00158836*freq)+ 1i*(0.00202237*freq)*((log(1/GMRii) + 7.6786 + 0.5*log(row/freq)));
                            else
                                 z(ii,jj,ConfgD(iii,1)) = 0  + (0.00158836*freq)+ 1i*(0.00202237*freq)*((log(1/Dis(ii,jj, ConfgD(iii,5))) + 7.6786 + 0.5*log(row/freq)));  
                            end
                        else
                            z(ii,jj,ConfgD(iii,1))=0;
                        end
                 end

            if  ConfgD(iii,10) == 1 % Concentric Neutral Cable  system   
              if (ConfgD(iii,jj+14)==1)&&(ConfgD(iii,ii+14)==1)  
                    if (ii==jj) && (jj<4)
                         [r2,c2,v2] = find(CondLibD(:,1)== ConfgD(iii,jj+5));% finding the index of the conductor lib in library data             
                         rii =  CondLibD(r2,4);
                         GMRii = CondLibD(r2,6);
                         z(ii,jj,ConfgD(iii,1)) = rii+(0.00158836*freq)+ 1i*(0.00202237*freq)*((log(1/GMRii) + 7.6786 + 0.5*log(row/freq)));

                    elseif (ii==jj) && (jj>=4)
                         [r2,c2,v2] = find(CondLibD(:,1)== ConfgD(iii,9));% finding the index of the new conductor in library data 
                         ds = CondLibD(r2,5); % diameter of neutral (valid for cable only!!)
                         rii =  CondLibD(r2,4)/k; % The equivalent resistance rcn = rs/k
                         R = (Dod-ds)/24; % Radius of the circle passing thru the cener of the strands
                         GMRii = CondLibD(r2,6);
                         GMRcn = (GMRii*k*R^(k-1))^(1/k);% Equivalent GMR of the concentric neutrals
                         z(ii,jj,ConfgD(iii,1)) = rii+(0.00158836*freq)+ 1i*(0.00202237*freq)*((log(1/GMRcn) + 7.6786 + 0.5*log(row/freq)));           

                    else 
                         z(ii,jj,ConfgD(iii,1)) = 0  + (0.00158836*freq)+ 1i*(0.00202237*freq)*((log(1/Dis(ii,jj,ConfgD(iii,5))) + 7.6786 + 0.5*log(row/freq)));  
                    end
              else
                  z(ii,jj,ConfgD(iii,1)) = 0;
              end
            end
        end
     end

     
% Tape shieled conductor
if  ConfgD(iii,10) == 2 % Tape Shielded Cable system      
       disp(' Tape shielded cable system is yet to be modelled')
       [r2,c2,v2] = find(CondLibD(:,1)== ConfgD(iii,6));
       [r2n,c2n,v2n] = find(CondLibD(:,1)==ConfgD(iii,9));

       ds =  ConfgD(iii, 12); % outside diameter
       rii = CondLibD(r2,4);
       GMRii = CondLibD(r2,6);
       T = ConfgD(iii, 13);
       rin = CondLibD(r2n,4);
       GMRin = CondLibD(r2n,6);
       Dis_pn = ConfgD(iii, 21);

       z(1,1,ConfgD(iii,1)) = 0.0953+rii+1i*(0.00202237*freq)*(log(1/GMRii)+7.6786 + 0.5*log(row/freq));
       z(1,2,ConfgD(iii,1)) = 0.0953+1i*(0.00202237*freq)*(log(1/(((ds/2)-(T/2000))/12))+7.6786 + 0.5*log(row/freq));  
       z(2,1,ConfgD(iii,1)) = z(1,2,ConfgD(iii,1));
       z(2,2,ConfgD(iii,1)) = 0.0953+(18.826/(ds*T))+1i*(0.00202237*freq)*(log(1/(((ds/2)-(T/2000))/12))+7.6786 + 0.5*log(row/freq));  
       z(1,3,ConfgD(iii,1)) = 0.0953+1i*(0.00202237*freq)*(log(1/(ConfgD(iii, 21)))+7.6786 + 0.5*log(row/freq)) ;
       z(3,1,ConfgD(iii,1)) = z(1,3,ConfgD(iii,1))  ;
       z(3,3,ConfgD(iii,1)) = 0.0953+rin+1i*(0.00202237*freq)*(log(1/GMRin)+7.6786 + 0.5*log(row/freq));
       z(2,3,ConfgD(iii,1)) = z(1,3,ConfgD(iii,1)) ;
       z(3,2,ConfgD(iii,1)) = z(1,3,ConfgD(iii,1)) ;

    B = z(1:1, 1:1,ConfgId(iii,1));
    C = z(1:1, 2:3,ConfgId(iii,1));
    D = z(2:3, 1:1,ConfgId(iii,1));
    E = z(2:3, 2:3,ConfgId(iii,1));
    z_abc(:,:,ConfgId(iii,1)) = B - C*inv(E)*D;
    z_abc(:,:,ConfgId(iii,1)) = z_abc(:,:,ConfgId(iii,1)).*[ConfgD(iii,15) 0 0; 0 ConfgD(iii,16) 0; 0 0 ConfgD(iii,17)];
end

if  ConfgD(iii,10) ~= 2
    B = z(1:3, 1:3,ConfgId(iii,1));
    C = z(1:3, 4:n,ConfgId(iii,1));
    D = z(4:n, 1:3,ConfgId(iii,1));
    E = z(4:n,4:n,ConfgId(iii,1));
    % Reduces the primitive matrices using the Kron Reduction technique
    z_abc(:,:,ConfgId(iii,1)) = B - C*inv(E)*D;
 end
     
end    
     
% End of Series Impedance Section
waitbar(0.2);
%% To find the shunt admittance of the line or cable

% Start of Shunt Admittance Section
y_abc = zeros(3,3);
for iii = 1:length(ConfgD(:,2))
    % To Compute the shunt admittance matrix

n = ConfgD(iii,11); % the primitive size of the matrix
Dod = ConfgD(iii,12);% Outside diameter (valid for cable only!!)
k = ConfgD(iii,13); % No of strands for cable type system 

if ConfgD(iii,10) == 0 % OH Transmission Line
  
  for ii=1:n
      for jj=1:n
        
[r2,c2,v2] = find(CondLibD(:,1)== ConfgD(iii,jj+5));% finding the index of the conductor lib in library data             
RD = CondLibD(r2,5)/24;
          if (ii==jj)
              P(ii,jj,ConfgId(iii,1)) = 11.17689*log(Sij(ii,jj,ConfgD(iii,5))/RD);             
          else
              P(ii,jj,ConfgId(iii,1)) = 11.17689*log(Sij(ii,jj,ConfgD(iii,5))/Dis(ii,jj,ConfgD(iii,5)));
          end
      end
  end
 
  B = P(1:3, 1:3,ConfgId(iii,1));
  C = P(1:3, 4:n,ConfgId(iii,1));
  D = P(4:n, 1:3,ConfgId(iii,1));
  E = P(4:n,4:n,ConfgId(iii,1));
      
 p_abc = B - C*inv(E)*D;
 
 % building the shunt admittance matrix
 %disp('The Shunt Capacitance Matrix is:')
C_abc = inv(p_abc);
 
 %The three phase shunt admittance matrix is given by Yabc = j*377*Cabc
 %disp('The Three Phase Shunt Admittance Matrix is:')
y_abc(:,:,ConfgId(iii,1)) = 1i*376.9911*C_abc;
  
end
    
  if ConfgD(iii,10) == 1 % Concentric Neutral Cable 
   [r3,c3,v3] = find(CondLibD(:,1) == ConfgD(iii,9)); % finding the index of the conductor lib in library data             
    RDs = CondLibD(r3,5)/2;
    ds = CondLibD(r3,5);
    R = (Dod-ds)/24; % Radius of the circle passing thru the cener of the strands
    Rb = R*12; % Radius in Inches
    %Y_abc = zeros(3,3);
        for ii=1:3
           [r2,c2,v2] = find(CondLibD(:,1) == ConfgD(iii,ii+5)); % finding the index of the conductor lib in library data             
            RDc = CondLibD(r2,5)/2;
            Yag = 0+1i*(77.3619/(log(Rb/RDc)-(1/k)*log(k*RDs/Rb)));
            y_abc(ii,ii,ConfgId(iii,1))  = Yag;
        end
  end
    
    
  for ii=1:3
    for jj=1:3
        if or(ConfgD(iii,ii+14)==0, ConfgD(iii,jj+14)==0)            
            y_abc(ii,jj,ConfgId(iii,1))=0;           
        end         
    end
  end
 
 
 
   if ConfgD(iii,10) == 2 % Tape shielded Cable 
   disp(' Need to model this section!');
   [r2,c2,v2] = find(CondLibD(:,1)== ConfgD(iii,6));
      
   ds =  ConfgD(iii, 12); % outside diameter
   Rdc = CondLibD(r2,5)/2;
   T = ConfgD(iii, 13);
   
   Rb = (ds-T/1000)/2;
   
   Yag = 0+1i*(77.3619/(log(Rb/Rdc)));
   y_abc(:,:,ConfgId(iii,1)) = Yag.*[ConfgD(iii,15) 0 0; 0 ConfgD(iii,16) 0; 0 0 ConfgD(iii,17)];

   end
 
end
disp(' End of Shunt Admittance Section')
waitbar(0.4);

%% UDB (User Defined Block creates the abcd parameter based on the Zabc & Yabc fed by user

for ii=1:length(UdbD(:,1))
z_abc(:,:,UdbD(ii,1)) = [UdbD(ii,2)+1i*UdbD(ii,3) UdbD(ii,4)+1i*UdbD(ii,5) UdbD(ii,6)+1i*UdbD(ii,7); 
                         UdbD(ii,4)+1i*UdbD(ii,5) UdbD(ii,8)+1i*UdbD(ii,9) UdbD(ii,10)+1i*UdbD(ii,11); 
                         UdbD(ii,6)+1i*UdbD(ii,7) UdbD(ii,10)+1i*UdbD(ii,11) UdbD(ii,12)+1i*UdbD(ii,13)];
                     
y_abc(:,:,UdbD(ii,1)) = 1i*[UdbD(ii,14) UdbD(ii,15) UdbD(ii,16); 
                            UdbD(ii,15) UdbD(ii,17) UdbD(ii,18); 
                            UdbD(ii,16) UdbD(ii,18) UdbD(ii,19)];
end

%% To generate the abcd and AB parameters for all the lines

Tableabcd = []
disp(' Generating the abcd parameters for Lines')
for ii=1:length(BrchD(:,2))
 % Follow the generalized abcd matrices 
U = diag(ones(3,1));

Z_abc = z_abc(:,:,BrchD(ii,5)).*(BrchD(ii,4)/5280);
Y_abc = y_abc(:,:,BrchD(ii,5)).*(BrchD(ii,4)*1e-6/5280);
a(:,:,BrchD(ii,1)) = U + 0.5*Z_abc*Y_abc;
b(:,:,BrchD(ii,1)) = Z_abc;
c(:,:,BrchD(ii,1)) = zeros(3,3);
d(:,:,BrchD(ii,1)) = U;
A(:,:,BrchD(ii,1)) = U;
B(:,:,BrchD(ii,1)) = inv(a(:,:,BrchD(ii,1)))*b(:,:,BrchD(ii,1));
disp('The generalized matrices are')

[r2,c2,v2] = find(ConfgD(:,1) == BrchD(ii,5));
for kk =1:3
    if ConfgD(r2,kk+14)==0
         a(kk,kk,BrchD(ii,1))=0;
         d(kk,kk,BrchD(ii,1))=0;
         A(kk,kk,BrchD(ii,1))=0;
         B(kk,kk,BrchD(ii,1))=0;    
    end 
end   

l_a(:,:,BrchD(ii,1)) = a(:,:,BrchD(ii,1)); 
l_b(:,:,BrchD(ii,1)) = b(:,:,BrchD(ii,1));
l_c(:,:,BrchD(ii,1)) = c(:,:,BrchD(ii,1)) ;
l_d(:,:,BrchD(ii,1)) = d(:,:,BrchD(ii,1)) ;
l_A(:,:,BrchD(ii,1)) = A(:,:,BrchD(ii,1));
l_B(:,:,BrchD(ii,1)) = B(:,:,BrchD(ii,1)) ;


% Added upon recommendation from Dr. Paudyal.
Tableabcd = [Tableabcd; 
    real(Z_abc)             , imag(Z_abc),... 
                              imag(Y_abc),...             
    real(l_a(:,:,BrchD(ii,1))), imag(l_a(:,:,BrchD(ii,1))),...
    real(l_b(:,:,BrchD(ii,1))), imag(l_b(:,:,BrchD(ii,1))),...
    real(l_c(:,:,BrchD(ii,1))), imag(l_c(:,:,BrchD(ii,1))),...
    real(l_d(:,:,BrchD(ii,1))), imag(l_d(:,:,BrchD(ii,1))),...
    real(l_A(:,:,BrchD(ii,1))), imag(l_A(:,:,BrchD(ii,1))),...
    real(l_B(:,:,BrchD(ii,1))), imag(l_B(:,:,BrchD(ii,1)))];

end
Brchs.Line.a = l_a;
Brchs.Line.b = l_b;
Brchs.Line.c = l_c;
Brchs.Line.d = l_d;
Brchs.Line.A = l_A;
Brchs.Line.B = l_B;

disp(' End of line abcd parameter estimation')

%% Generation of abcd parameters for the Transformer models
 % BrchNo--> No of Branches , BrchId --> Branch Ids
 disp(' Start of Trf abcd parameter estimation');

D = [1 -1  0; 
     0  1 -1; 
    -1  0  1];

disp(' Start of Trf abcd parameter estimation');
TrfNo = length(find(TrfD(:,1)));
if TrfNo ~= 0 % if no trf data is available then skip this code

TrfId = TrfD(:,1);
TrfKva = TrfD(:,5);

 for iii = 1:TrfNo % for checking the complete trf datas
    if (TrfD(iii,7)==0 && TrfD(iii,9)~= 0) % Delta-Wye grounded/ungrounded
        TrfNt(iii,1) = TrfD(iii,6)/(TrfD(iii,8)/sqrt(3)); % nt = kVLL_pri/kVLN_sec 

    elseif (TrfD(iii,7)==1 && TrfD(iii,9)==0) % Ungrounded Wye-Delta 
        TrfNt(iii,1) = (TrfD(iii,6)/sqrt(3))/TrfD(iii,8);% nt = kVLN_pri/kVLL_sec        

    elseif (TrfD(iii,7)==2 && TrfD(iii,9)==0) % Grounded Wye-Delta
        TrfNt(iii,1) = (TrfD(iii,6)/sqrt(3))/TrfD(iii,8);% nt = kVLN_pri/kVLL_sec

    elseif (TrfD(iii,7)==0 && TrfD(iii,9)==0) % Delta-Delta
        TrfNt(iii,1) = (TrfD(iii,6)/TrfD(iii,8));% nt = kVLL_pri/kVLL_sec

    else
        TrfNt(iii,1) = (TrfD(iii,6)/TrfD(iii,8));
    end
       

% <<<<<<<<< Calculation of Zbase and Transformer impedances >>>>>>>>> 
    if TrfD(iii,13)== 0 % Trf kVA as total kVA for all phases  
        % Calculating the base impedances
        if TrfD(iii,9)==0 % Delta configuration on LV side
        TrfZbase(iii,1) = TrfD(iii,8)*TrfD(iii,8)*1000/(TrfD(iii,5)/3);    
            
        elseif TrfD(iii,9)==2 % WYE grounded on LV side 
        TrfZbase(iii,1) = (TrfD(iii,8)/sqrt(3)).^2*1000/(TrfD(iii,5)/3); 
        
        else % WYE ungrounded on LV side
        TrfZbase(iii,1) = TrfD(iii,8).^2*1000/(TrfD(iii,5)/3);    
        end
       
        % The Transformer impedance referenced to the low voltage side is
        Zt(iii,1) = ((TrfD(iii,10)+i*TrfD(iii,11))/100)*TrfZbase(iii,1);
        Ztabc(:,:,iii) = Zt(iii,1).*eye(3,3);% The transformer phase impedance matrix is
            
     % Trf kVA as total kVA for all phases     
    elseif TrfD(iii,13)==1 % kVA defined as separately for each phases
        TrfZbase(iii,1) = TrfD(iii,8)*TrfD(iii,8) * 1000/TrfD(iii,14); 
        TrfZbase(iii,2) = TrfD(iii,8)*TrfD(iii,8) * 1000/TrfD(iii,15); 
        TrfZbase(iii,3) = TrfD(iii,8)*TrfD(iii,8) * 1000/TrfD(iii,16); 

      % the Transformer impedance referenced to the low voltage side is
        Zt(iii,1) = ((TrfD(iii,17)+i*TrfD(iii,18))/100)*TrfZbase(iii,1);  
        Zt(iii,2) = ((TrfD(iii,19)+i*TrfD(iii,20))/100)*TrfZbase(iii,2);
        Zt(iii,3) = ((TrfD(iii,21)+i*TrfD(iii,22))/100)*TrfZbase(iii,3);

      % The transformer phase impedance matrix is
        Ztabc(:,:,iii) = [Zt(iii,1) 0 0 ; 0 Zt(iii,2) 0; 0 0 Zt(iii,3)];

    elseif TrfD(iii,13)== 3 % Impedance computed by user 
    Ztabc(:,:,iii) = [UdbT(iii,2)+i*UdbT(iii,3) 0 0 ; 
                       0 UdbT(iii,8)+i*UdbT(iii,9) 0; 
                       0 0 UdbT(iii,12)+i*UdbT(iii,13)];
    end
    
    
% <<<<<<<<<<<<<<<< Computation of abcd parameters >>>>>>>>>>>>>>>>>>>>   
% Delta - Grounded Wye connection
if (TrfD(iii,7)==0 && TrfD(iii,9)==2)  
    disp('Delta - Grounded Wye connection')
    t_a(:,:,TrfD(iii,1)) = -TrfNt(iii,1)/3 .*[ 0 2 1; 1 0 2; 2 1 0];
    t_b(:,:,TrfD(iii,1)) = t_a(:,:,TrfD(iii,1))*Ztabc(:,:,iii);
    t_A(:,:,TrfD(iii,1)) = 1/TrfNt(iii,1) .* [1 0 -1; -1 1 0; 0 -1 1];
    t_B(:,:,TrfD(iii,1)) = Ztabc(:,:,iii);
    t_c(:,:,TrfD(iii,1)) = [0 0 0;0 0 0;0 0 0];
    t_d(:,:,TrfD(iii,1)) = 1/TrfNt(iii,1) .* [1 -1 0; 0 1 -1; -1 0 1];        
end

% Grounded Wye-Grounded Wye Step down connection   
if (TrfD(iii,7)==2 && TrfD(iii,9)==2)   
    disp('Grounded Wye-Grounded Wye Step down connection')
    t_a(:,:,TrfD(iii,1)) = TrfNt(iii,1).*[ 1 0 0; 0 1 0; 0 0 1];
    t_b(:,:,TrfD(iii,1)) = t_a(:,:,TrfD(iii,1))*Ztabc(:,:,iii);
    t_A(:,:,TrfD(iii,1)) = inv(t_a(:,:,TrfD(iii,1)));
    t_B(:,:,TrfD(iii,1)) = Ztabc(:,:,iii);
    t_c(:,:,TrfD(iii,1)) = [0 0 0;0 0 0;0 0 0];
    t_d(:,:,TrfD(iii,1)) = 1/TrfNt(iii,1) .* [1 0 0; 0 1 0; 0 0 1];
end

%Delta-Delta Connection
if (TrfD(iii,7)==0 && TrfD(iii,9)==0)
    disp('Delta-Delta Connection')        
    F =[ 1                  0               -1; 
        -1                  1                0;
        Ztabc(1,1,iii) Ztabc(2,2,iii) Ztabc(3,3,iii)];
    G = inv(F);  G(1:3,3:3)=0;  G1 = G;
    AV = TrfNt(iii,1)*[ 1 0 0; 
                        0 1 0; 
                        0 0 1];                    
    AI = TrfNt(iii,1)*[ 1 0 0; 
                        0 1 0; 
                        0 0 1];  

    t_A(:,:,TrfD(iii,1)) = (1/(3*TrfNt(iii,1))).*[2 -1 -1;-1 2 -1;-1 -1 2] ;
    t_B(:,:,TrfD(iii,1)) = W*Ztabc(:,:,iii)*G1;
    t_a(:,:,TrfD(iii,1)) = W*AV*D;
    t_b(:,:,TrfD(iii,1)) = W*AV*Ztabc(:,:,iii)*G1;
    t_c(:,:,TrfD(iii,1)) = [0 0 0;0 0 0;0 0 0];
    t_d(:,:,TrfD(iii,1)) = 1/TrfNt(iii,1) .* [1 0 0; 0 1 0; 0 0 1];
end

%Grounded Wye and Delta Step down connection
if (TrfD(iii,7)==2 && TrfD(iii,9)==0) 
    disp('Grounded Wye and Delta Step down connection')
      
F =[ 1                  0                -1; 
    -1                  1                 0;
    TrfNt(iii,1)*Ztabc(1,1,iii)+((3/TrfNt(iii,1))*TrfD(iii,23))...
    TrfNt(iii,1)*Ztabc(2,2,iii)+((3/TrfNt(iii,1))*TrfD(iii,23)) ...
    TrfNt(iii,1)*Ztabc(3,3,iii)+((3/TrfNt(iii,1))*TrfD(iii,23))];

    G = inv(F);   G1 = [G(:,3) G(:,3) G(:,3) ];  G2=G;   G2(:,3)=0;
    
    AV = TrfNt(iii,1)*[ 1 0 0; 
                        0 1 0; 
                        0 0 1];
                    
    AI = TrfNt(iii,1)*[ 1 0 0; 
                        0 1 0; 
                        0 0 1];
                    
    Zg = TrfD(iii,23)*ones(3,3);                
    X1 = Ztabc*AI + inv(AV)*Zg;                   

    t_a(:,:,TrfD(iii,1)) = zeros(3,3);
    t_b(:,:,TrfD(iii,1)) =  zeros(3,3);    
    t_c(:,:,TrfD(iii,1)) = inv(AI)*G1;
    t_d(:,:,TrfD(iii,1)) = inv(AI)*G2;
    t_A(:,:,TrfD(iii,1)) = W*(inv(AV)-X1*t_c(:,:, TrfD(iii,1)));
    t_B(:,:,TrfD(iii,1)) = W*X1*t_d(:,:,TrfD(iii,1));    
end

if (TrfD(iii,7)==1 && TrfD(iii,9)==0)    % Ungrounded Wye-Delta Step down connection
     disp('Ungrounded Wye-Delta Step down connection');    
    t_a(:,:,TrfD(iii,1)) = TrfNt(iii,1) .*[ 1 -1  0; 
                                            0  1 -1; 
                                           -1  0  1]      
    t_b(:,:,TrfD(iii,1)) = TrfNt(iii,1)/3 .*...
                  [Ztabc(1,1,iii), -Ztabc(1,1,iii),   0;
                  Ztabc(2,2,iii),   2*Ztabc(2,2,iii), 0;
                 -2*Ztabc(3,3,iii), -Ztabc(3,3,iii),  0]
    t_A(:,:,TrfD(iii,1)) = 1/(3*TrfNt(iii,1)) .* [2 1 0; 
                                                  0 2 1; 
                                                  1 0 2]
    t_B(:,:,TrfD(iii,1)) = 1/9.*...
    [2*Ztabc(1,1,iii)+Ztabc(2,2,iii), 2*Ztabc(2,2,iii)-2*Ztabc(1,1,iii), 0;
     2*Ztabc(2,2,iii)-2*Ztabc(3,3,iii), 4*Ztabc(2,2,iii)-Ztabc(3,3,iii), 0;
     Ztabc(1,1,iii)-4*Ztabc(3,3,iii), -Ztabc(1,1,iii)-2*Ztabc(3,3,iii), 0]
    
    t_c(:,:,TrfD(iii,1)) = [0 0 0;
                            0 0 0;
                            0 0 0];
    t_d(:,:,TrfD(iii,1)) = 1/(3*TrfNt(iii,1)) .* [1 -1 0; 
                                                  1  2 0; 
                                                 -2 -1 0]
 end
 
 Tableabcd = [Tableabcd;real(Ztabc(:,:,iii)), imag(Ztabc(:,:,iii)),...
     zeros(3,3),...
     real(t_a(:,:,TrfD(iii,1))), imag(a(:,:,TrfD(iii,1))),...
     real(t_b(:,:,TrfD(iii,1))), imag(t_b(:,:,TrfD(iii,1))),...
     real(t_c(:,:,TrfD(iii,1))), imag(t_c(:,:,TrfD(iii,1))),...
     real(t_d(:,:,TrfD(iii,1))), imag(t_d(:,:,TrfD(iii,1))),...
     real(t_A(:,:,TrfD(iii,1))), imag(t_A(:,:,TrfD(iii,1))),...
     real(t_B(:,:,TrfD(iii,1))), imag(t_B(:,:,TrfD(iii,1)))];
 
   end
else 
    disp('Transformer data is not available'); 
    t_a = [];
    t_b = [];
    t_c = [];
    t_d = [];
    t_A = [];
    t_B = [];
    TrfNt = [];
    TrfNo = [];
    Ztabc = [];
end

Brchs.Trf.a = t_a;
Brchs.Trf.b = t_b;
Brchs.Trf.c = t_c;
Brchs.Trf.d = t_d;
Brchs.Trf.A = t_A;
Brchs.Trf.B = t_B;

% TBrchs

disp(' End of trf abcd parameter estimation');
disp('end');

%% Generation of abcd parameters for the Voltage Regulator Section
if length(find(RegD(:,1)))~=0 && RegD(1,1)~=0    
%     TapInfo = [10;8;11]
iter=0
    [Tap, Brchs] = VoltReg(RegD, BusD, Brchs, iter);
end

%%
clk = clock;
clk1 = fix(clk)
disp('End of Input file generation');
disp(' <<                Process Completed                           >>')

waitbar(1);
close(stp1);
length(CapsD);
%%

% save('Tableabcd.mat','Tableabcd')
save ('inputdata.mat','Brchs','input')

save InputFileGenerator.mat ;
% this section is used to place all the .mat files in resp folder
% Binfile=fullfile(recycle, 'inputdata.mat');
% save (Binfile);
% 
% Binfile=fullfile(recycle, 'command.mat');
% save (Binfile);



