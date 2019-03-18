function [v,I, IL, ILpq, ILc, ILz, info] = ThreePhLF(hr, SCase, CapsD, RegD)
%% TreeTab function execution
   %run Input_file_generatorR8.m
  
% 1 - Snapshot Power flow (static load data type);
% 0- Time-series load data   
loadflag = 1;  
load Inputdata.mat ;
load command.mat ;
% User input for process
% SCase = 1;


load('TreeTab.mat')
Tree = array2table(TreeTab);
%     Tree.Properties.VariableNames = {'Si_no','BusSend','BrchType','BusRecv','c5','ElemID1','c7','ElemID2','c9'} 
    
    %% Three Phase load flow

%load InputFileGenerator.mat
tes=tic;
load TreeTab.mat ;
%load LFS_data.mat  ;
DI = [1 0 -1; -1 1 0; 0 -1 1];
D = [1 -1 0; 0 1 -1; -1 0 1];
W = 1/3*[2 1 0; 0 2 1; 1 0 2];




%D=[1 1 1;1 1 1;1 1 1];
Bc=zeros(3,1);
tol = input.data.StudyCase(SCase,3);
iter = 0; val = 99999; tol = 1e-4;
%% loading info
BusD = input.data.Nodes;
BrchD = input.data.Branch; % Branch Data imports
ConfgD = input.data.Configuration ;% Branch data configuration import
LoadsD = input.data.Loads; % Load data import

% if loadflag == 0
% % ****** Load data is available from IEEE and stored in *.mat format here.
% load Euro_LV_CSV.mat
% for kl =1:55
%     LoadsD(kl,4) = data(hr(1),2,kl) ;
%     LoadsD(kl,6) = data(hr(1),2,kl) ; 
%     LoadsD(kl,8) = data(hr(1),2,kl) ;
% end
% end


TrfD = input.data.Transformers; % Transformer data import
% RegD = input.data.Regulator; % Regulator data import
StudyInfo = input.data.StudyCase; % STudy cases imports


%% Initializing the matrices for the analysis
% Initializing all nodes to voltage matrix

nb = length(find(BusD(:,1)));
for ii = 1:nb
v(1:3,1:1,BusD(ii,1)) = 0+0*1i;
V_old(1:3,1:1,BusD(ii,1)) = 0+0*1i;
end

% Initializing the source voltage matrix
% src = [701 702]
n = length(src); % No of source feeder
for iisrc = 1:length(src)


%  v(:,:,src(1,iisrc)) = [conver2rec(6665.09, -0.1);
%     conver2rec(6662.28,-120);
%     conver2rec(6667.25, 120)];


% v(:,:,src(1,iisrc)) = [conver2rec(12.66e3, 0);
%     conver2rec(12.66e3, -120);
%     conver2rec(12.66e3, 120)];

%  v(:,:,src(1,iisrc)) = [conver2rec(251.901, -30.2);
%     conver2rec(251.443, -150.3);
%     conver2rec(251.952, 89.9)];
    


% % k=[conver2rec(1.021,   -2.49 );conver2rec(1.042,   -1.72 ); conver2rec(1.017,  -2.17 )]
k=[conver2rec(1.00,   0 );conver2rec(1.00,   0 ); conver2rec(1.00,  0 )]

[isrc,~,~] = find(BusD(:,1) == src(1,iisrc));


% % ******   IEEE 13 bus system; Voltage at Node 650  ************
 v(:,:,src(1,iisrc)) =   [conver2rec(BusD(isrc,4)*577.3503,   0.000 );
                       conver2rec(BusD(isrc,4)*577.3503,  -120.00);
                       conver2rec(BusD(isrc,4)*577.3503, 120.00)]; 


% % % ******   IEEE 13 bus system; Voltage at Node 650  ************
%  v(:,:,src(1,iisrc)) =   [conver2rec(2401.777,   0.000 );
%                        conver2rec(2401.777,  -120.00);
%                        conver2rec(2401.777, 120.00)];   



% % ****************** Single Phase 730 Node System 
% k=[conver2rec(1.00,   0 );conver2rec(1.00,   0 ); conver2rec(1.00,  0 )];
% 
% v(:,:,src(1,iisrc)) = [conver2rec(12.66e3, 0);
%     conver2rec(12.66e3, -120);
%     conver2rec(12.66e3, 120)];                   
%                    


% % *****************  IEEE European system Time series 
% k=[conver2rec(1.05,   0 );conver2rec(1.05,   0 ); conver2rec(1.05,  0 )];
% 
%  v(:,:,src(1,iisrc)) = W*[conver2rec(11e3, 30);
%     conver2rec(11e3, -90);
%     conver2rec(11e3, 150)];                   
                   
                   
end

 % Intializing the current matrices for Lines

% n = length(find(BrchD(:,1)));
% for ii = 1:n  
% I(:,:,BrchD(ii,1),1) = zeros(3,1); %0+0*1i;                  
% end

n = max(BrchD(:,1));
if n == 0
    n = 1;
end
I(:,:,n,1) = zeros(3,1);  % Intializing the current matrices for Lines


n = max(TrfD(:,1));
if n == 0
    n = 1;
end

I(:,:,n,2) = zeros(3,1);  % Intializing the current matrices for Trf


n=max(RegD(:,1));
if n==0
    n=1;
end
I(:,:,n,3) = zeros(3,1);  % Intializing the current matrices for Regulator


% Intializing the current matrices for Load
% Intializing the current matrices for Bus (For summing the shunt element
% current)
% Intializing the current matrices for capacitor
% n = length(find(BusD(:,1)));
for ii = 1:length(find(BusD(:,1)))  
 IL(1:3,1:1,BusD(ii,1)) = 0;  
 Ibus(1:3,1:1,BusD(ii,1)) = 0; 
 IC(1:3,1:1,BusD(ii,1)) = 0;
 I_abc(1:3,1:1,BusD(ii,1))=0;
 ILall(1:3,1:3:BusD(ii,1))=0;
end

%n = length(find(LoadsD(:,1)));
for ii = 1:length(find(LoadsD(:,1)))  
 I(1:3,1:1,LoadsD(ii,1),8) = 0;  
end

% n = length(find(CapsD(:,1)));
for ii = 1:length(find(CapsD(:,1)))  
 I(1:3,1:1,CapsD(ii,1),9) = 0;  
end


% Regulator data processing
if RegD(1,1)~=0
[Tap, Brchs] = VoltReg(RegD, BusD, Brchs, iter);
end



% load Euro_LV_CSV.mat

%% Backward Sweep Technique
for h=1:length(hr)

    tBFS = tic;  
% if loadflag == 0
% % ****** Load data is available from IEEE and stored in *.mat format here.
% 
% for kl =1:55
%     LoadsD(kl,4) = data(hr,2,kl) ;
%     LoadsD(kl,6) = data(hr,2,kl) ; 
%     LoadsD(kl,8) = data(hr,2,kl) ;
% end
% end


for ii = 1:length(find(BusD(:,1)))  
 IL(1:3,1:1,BusD(ii,1)) = 0;  
 Ibus(1:3,1:1,BusD(ii,1)) = 0; 
 IC(1:3,1:1,BusD(ii,1)) = 0;
 I_abc(1:3,1:1,BusD(ii,1))=0;
 ILall(1:3,1:3:BusD(ii,1))=0;
end

%n = length(find(LoadsD(:,1)));
for ii = 1:length(find(LoadsD(:,1)))  
 I(1:3,1:1,LoadsD(ii,1),8) = 0;  
end

% n = length(find(CapsD(:,1)));
for ii = 1:length(find(CapsD(:,1)))  
 I(1:3,1:1,CapsD(ii,1),9) = 0;  
end


iter = 0;
 val=999999;
 
 
 
while (tol<val) && iter<500  %StudyInfo(SCase,4)
%     (iter~=0 & tol<val) && iter<500
iter  ;
m =  length(TreeTab(:,1)); %No of rows in a TreeTable Matrix


% creating time loop for time-series load flow


    
for ii = 1:m
       nii = length(find(TreeTab(ii,4:3+BrchMax)));% counts the no of branches for a node that needs execution
        for jj = 1:nii            
            % Associate the proper A and B Matrices of the
            % line/trf/regulator            
            if (TreeTab(ii,BrchMax+3+jj)== 1)% For lines/cables
                [r3,c3,v3] = find(BrchD(:,1)==TreeTab(ii,2*BrchMax+3+jj));
                if (BrchD(r3,6)==1)
                A = Brchs.Line.A(:,:,TreeTab(ii,2*BrchMax+3+jj));
                B = Brchs.Line.B(:,:,TreeTab(ii,2*BrchMax+3+jj));                
                %A = l_A(:,:,TreeTab(ii,2*BrchMax+3+jj));
                %B = l_B(:,:,TreeTab(ii,2*BrchMax+3+jj));                
                I_ABC = I(:,:,TreeTab(ii,2*BrchMax+3+jj),1);
%                 I_ABC = Ibus(:,:,TreeTab(ii,3+jj))
                end
            end            
            if TreeTab(ii,BrchMax+3+jj)== 2 % For Trfs
%                 A = t_A(:,:,TreeTab(ii,2*BrchMax+3+jj));
%                 B = t_B(:,:,TreeTab(ii,2*BrchMax+3+jj));                
                A = Brchs.Trf.A(:,:,TreeTab(ii,2*BrchMax+3+jj));
                B = Brchs.Trf.B(:,:,TreeTab(ii,2*BrchMax+3+jj))  ;            
                               
%                  I_ABC = I(:,:,TreeTab(ii,2*BrchMax+3+jj),2)
                I_ABC = Ibus(:,:,TreeTab(ii,3+jj));
            end            
            if TreeTab(ii,BrchMax+3+jj)==3 % )%&&(RegD(TreeTab(ii,2*BrchMax+3+jj),29)==0)  % for voltage regulator
%                [Tap, Brchs] = VoltReg(RegD,BusD,iter,v,Ibus,[],I);
                if iter==0
                    A = eye(3,3);
                else
                    A = Brchs.Reg.A(:,:,TreeTab(ii,2*BrchMax+3+jj));
                end
                B = Brchs.Reg.B(:,:,TreeTab(ii,2*BrchMax+3+jj));
                I_ABC = I(:,:,TreeTab(ii,2*BrchMax+3+jj),3);
            end
            v(:,:,TreeTab(ii,3+jj)) = A*v(:,:,TreeTab(ii,2))- B*I_ABC;           
        end         
end 
 
 %%
 % |*Convergence checker & compute the load current and Shunt capacitor.*|
 % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   
  for ii = 1:length(find(TermBus(:,1)))
  [r2,~,~] = find(BusD(:,1) == TermBus(ii,1)) ;
  Er(1:3,1:1,TermBus(ii,1)) = abs(abs(v(:,:,TermBus(ii,1)))-abs(V_old(:,:,TermBus(ii,1))))./(BusD(r2,4)*1000/sqrt(3)) ;
  mz(ii,1) = max(Er(:,:,TermBus(ii,1)));
       
  end
 
  % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
  
 val = max(mz(:,1));
  
 % Compute Load current 
 % Adding all the phases of the load

% [Ibus,ILpq,ILz,ILc,IC] = LoadnCap(v,iter,IC,IL, LoadsD);



% if loadflag == 0
% % ****** Load data is available from IEEE and stored in *.mat format here.
% load Euro_LV_CSV.mat
% for kl =1:55
%     LoadsD(kl,4) = data(hr(h),2,kl) ;
%     LoadsD(kl,6) = data(hr(h),2,kl) ; 
%     LoadsD(kl,8) = data(hr(h),2,kl) ;
% end
% end



[Ibus,IL, ILpq, ILz, ILc] = LoadnCap(v,iter,IL, LoadsD); % Load current calc
    if iter==1
        n = length(src); % No of source feeder
        for iisrc = 1:length(n)
       v(:,:,src(1,iisrc)) =  k.*v(:,:,src(1,iisrc));
        end
    end
[Ibus, IC] = Caps(v,IC,IL,Ibus, CapsD,BusD); % Capacitor current calc

iter = iter+1;

% clear I_abc Matrix before KCL else might result in dubious results

for ii = 1:length(BusD(:,1)) 
 I_abc(1:3,1:1,BusD(ii,1))=0;
end

%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 %                       
  % *Forward Sweep Technique*
 %
 % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 m = length(TreeTab(:,1));

 for ii = m:-1:1     
    for jj = 1:TreeTab(ii,3) % looks for no of branches and executes 
        % each branch separately thereby reducing the time consumed 
        if TreeTab(ii,3+jj)~=0 % only non terminal bus and branches would be executed.

            if (TreeTab(ii,3+BrchMax+jj) == 1)%&&(BrchD(TreeTab(ii,2*BrchMax+3+jj),6)==1) % if the branch is a Line & if inservice only
                
                d = Brchs.Line.d(:,:,TreeTab(ii,3+(2*BrchMax)+jj));
                
                I_abc(:,:,TreeTab(ii,2)) = d*Ibus(:,:,TreeTab(ii,3+jj))+I_abc(:,:,TreeTab(ii,2));
                I(:,:,TreeTab(ii,2*BrchMax+3+jj),1) = d*Ibus(:,:,TreeTab(ii,3+jj));
            end
                    
                    if TreeTab(ii,3+BrchMax+jj) == 2 % if the branch is a transformer
                        
                        d = Brchs.Trf.d(:,:,TreeTab(ii,3+(2*BrchMax)+jj));
                        
                        I_abc(:,:,TreeTab(ii,2)) = d*Ibus(:,:,TreeTab(ii,3+jj))+I_abc(:,:,TreeTab(ii,2));
                        I(:,:,TreeTab(ii,2*BrchMax+3+jj),2) = d*Ibus(:,:,TreeTab(ii,3+jj));
                    end
                    
            if (TreeTab(ii,3+BrchMax+jj) == 3) %&&(RegD(TreeTab(ii,3+(2*BrchMax)+jj),29)==1) % if the branch is a Regulator
%                 [Tap,r_a,r_d,r_A,r_B] = VoltReg(RegD,BusD,iter,v,Ibus,[],I);
                d = Brchs.Reg.d(:,:,TreeTab(ii,3+(2*BrchMax)+jj));
%                 d = r_d(:,:,TreeTab(ii,3+(2*BrchMax)+jj));
                
                I_abc(:,:,TreeTab(ii,2)) = d*Ibus(:,:,TreeTab(ii,3+jj))+I_abc(:,:,TreeTab(ii,2));
                I(:,:,TreeTab(ii,2*BrchMax+3+jj),3) = d*Ibus(:,:,TreeTab(ii,3+jj));
            end
                    
    % Note: Here the series elements are at atmost imporatant 
          % since all the parallel/shunt element current are computed and 
          % stored in Ibus matrix                    
        end
          
%Ibus(:,:,TreeTab(ii,2))=Ibus(:,:,TreeTab(ii,2))+I_abc(:,:,TreeTab(ii,2));
% KCL at the sending end bus is carried out here

    end

  
  if TreeTab(ii,2)~=0
    Ibus(:,:,TreeTab(ii,2))=Ibus(:,:,TreeTab(ii,2))+I_abc(:,:,TreeTab(ii,2));  
  end
  
 end
 
      for ii = 1:length(BusD(:,1))
        V_old(:,:,BusD(ii,1)) = v(:,:,BusD(ii,1));
      end
     
end


% Vnode(:,:,h)= v(:,1,34);
% % abs(Vnode);
% Inode(:,:,:, h)=Ibus;
% %Iteration counts;
% PowerFlowIter (h,1) = iter;
% % disp('End of code')     
time(h,1) = toc(tBFS);

% % m =  length(TreeTab(:,1)); 
% % for ii = 1:length(TreeTab(:,1)) %No of rows in a TreeTable Matrix
% % 
% %        nii = length(find(TreeTab(ii,4:3+BrchMax)));% counts the no of branches for a node that needs execution
% %         for jj = 1:TreeTab(ii,3)
% %            if TreeTab(ii,3+jj+BrchMax)~=0;
% %                
% %             if TreeTab(ii,BrchMax+3+jj)== 1 % For Lines and Cables
% %                 
% %                 A = Brchs.Line.A(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 B = Brchs.Line.B(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 C = Brchs.Line.c(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 D = Brchs.Line.d(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 
% %             end
% %                        
% %             if TreeTab(ii,BrchMax+3+jj)== 2 % For Trfs
% %                 
% %                 A = Brchs.Trf.A(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 B = Brchs.Trf.B(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 C = Brchs.Trf.c(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 D = Brchs.Trf.d(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %             
% %             end
% %             
% %             if TreeTab(ii,BrchMax+3+jj)==3 % For Regulator
% %                 A = Brchs.Reg.A(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 B = Brchs.Reg.B(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 C = Brchs.Reg.c(:,:,TreeTab(ii,2*BrchMax+3+jj));
% %                 D = Brchs.Reg.d(:,:,TreeTab(ii,2*BrchMax+3+jj)) ;               
% %             end
% %             
% % % Receiving end current
% % Ibrch(:,2,TreeTab(ii,BrchMax+3+jj),TreeTab(ii,2*BrchMax+3+jj))= inv(B)*(v(:,:,TreeTab(ii,2))-A*v(:,:,TreeTab(ii,3+jj)));
% % 
% % % Sending end current
% % Ibrch(:,1,TreeTab(ii,BrchMax+3+jj),TreeTab(ii,2*BrchMax+3+jj))= C*v(:,:,TreeTab(ii,3+jj))+D*Ibrch(:,2,TreeTab(ii,BrchMax+3+jj),TreeTab(ii,2*BrchMax+3+jj));       
% %            end            
% %         end
% % end

end
% disp(time)

info.exetime = time;
info.iter=iter;

%%
clear A B a b c d C D Bc IC ii jj m n nb nii r2 c2 v2 r3 %Brchs 
clear tser 
clear functions
save LFSresult.mat
% disp(abs(v(1,1,34)))
% disp(LoadsD(1,4))
% disp(time)
% disp(iter)

% disp (' ****************** END ************************')

end

% TreeTab
% 
%%To form the tree table
%TreeTab = [SiNo BusNo,BrchNo, RcvBus1, RcvBus2, ...... ,RcvBusN, Typ1, Typ2..
%....,TypN, ElemId1, ElemId2, ElemId2,......ElemIdN] 
