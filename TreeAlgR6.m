function [TreeTab] = TreeAlgR6(SCase)
% load InputFileGenerator.mat
TermBus=[];

load Inputdata.mat;
BusD = input.data.Nodes;
BrchD = input.data.Branch; % Branch Data imports
ConfgD = input.data.Configuration ;% Branch data configuration import
LoadsD = input.data.Loads; % Load data import
TrfD = input.data.Transformers; % Transformer data import
RegD = input.data.Regulator; % Regulator data import
CapsD = input.data.Caps; % Shunt capacitor details import
StudyInfo = input.data.StudyCase; % STudy cases imports


TrfNo = length(find(TrfD(:,1)));
RegNo = length(find(RegD(:,1)));


    % src = input('Enter the source bus No:');
    % src is like a slack bus in a system.

    src = StudyInfo(SCase,9) % SOurce Bus
    n = length(BusD(:,1)); % no of buses
    nsrc = length(src); % no of sources
    co=0;
    for ii = 1:n
        if(BusD(ii,1)==src);
            co = co+1;
        end
    end
    if co == 0
        X= ['Bus No: ',num2str(src), '  is not available in the bus data, restart and enter',' the correct source bus'];
        disp(X)
    end
    BusNo = length(BusD(:,1));
    in = 0;

    % Count the No of branches for a given Node.
    % Counts the number of lines for a node by looking into from and to bus
    % of the branches
    for ii = 1:BusNo
        for jj = 1:length(BrchD(:,2))
        if BusD(ii,1) == BrchD(jj,2) %/ BusD(ii,1) == BrchD(jj,3)
            in = in+1;
        end
         if BusD(ii,1) == BrchD(jj,3) %/ BusD(ii,1) == BrchD(jj,3)
            in = in+1;
        end
        end
        BusChk(ii,1)=BusD(ii,1);
        BusChk(ii,2)=in;
        in = 0;
    end

    %disp('counts the no of transformer terminals at a given node')
    for ii = 1:BusNo
        for jj = 1:TrfNo
        if BusD(ii,1) == TrfD(jj,3) %/ BusD(ii,1) == BrchD(jj,3)
            in = in+1;
        end
         if BusD(ii,1) == TrfD(jj,4) %/ BusD(ii,1) == BrchD(jj,3)
            in = in+1;
        end
        end

        BusChk(ii,2)=  BusChk(ii,2)+in;
        in = 0;
    end


    %counts the no of switches at a given node


    % disp('counts the no of regulator terminals at a given node')
    for ii = 1:BusNo
        for jj = 1:RegNo
        if BusD(ii,1) == RegD(jj,2) %/ BusD(ii,1) == BrchD(jj,3)
            in = in+1;
        end
         if BusD(ii,1) == RegD(jj,3) %/ BusD(ii,1) == BrchD(jj,3)
            in = in+1;
        end
        end
       %BusChk(ii,3)=  in;
       BusChk(ii,2)=  BusChk(ii,2)+in;
        in = 0;
        BusChk(ii,3) = 0;
    end


    % To find the load bus i.e., terminal bus for a node.
    %disp('The following code generates only the node that has 1 branch and a load.')
    co=1;
    for ii = 1:length(LoadsD(:,1))
        for jj = 1:BusNo
            if LoadsD(ii,1)== BusD(jj,1)&& BusChk(jj,2)==1
                TermBus(co,1) = BusD(jj,1);
                co=co+1;
            end
        end
    end

    % placing 1 for slack bus (intializing)
    for ii = 1:BusNo
        if BusChk(ii,1) == src
        BusChk(ii,3) = 1;
        end
    end

    % Initilizes the Branch data if the from or to bus has the slack bus
    for ii = 1:length(BrchD(:,2))
        if BrchD(ii,2)==src
            BrchD(ii,6) = -1;
        end
            if BrchD(ii,3)==src
            BrchD(ii,7) = -1;
            end      
    end


    %%To form the tree table
    %TreeTab = [SiNo BusNo,BrchNo, RcvBus1, RcvBus2, ...... ,RcvBusN, Typ1, Typ2..
    %....,TypN, ElemId1, ElemId2, ElemId2,......ElemIdN] 
    disp(' Start of Tree Table')
    BrchMax = max(BusChk(:,2));
    Typ=BrchMax+3;
    Elem=2*BrchMax+3;
    TreeTab = zeros(BusNo,3*BrchMax+3);
    Temp=zeros(1,1);
    NoTemp = 0;
    co = 1; % Column index
    si = 1;

    ln = length(src);
    if ln > 1 
        for ii = 1:ln
            ptr = src(1,ii);
        end  

        for jj = 1:ln
            Temp(1,jj) = src(1,jj+1);
            NoTemp = jj;
        end
    end
    TreeTab(1,2)= src(1,1);

    disp('start of tree algor')
    status = 0;
    CurBus = TreeTab(si,2);
    %while(status == 1)
    for nnn = 1:10000
    N = 0;
    TreeTab(si,1) = si; % SI Column
    TreeTab(si,2) = CurBus ;

    [r3,c3,v3] = find(BusChk(:,1)==CurBus);
    TreeTab(si,3) = BusChk(r3,2);

    %for ii = 1:TreeTab(si,3)
        co = 4;
    % finding the nodes from the Lines
    [r,c,v] = find(BrchD(:,2)==CurBus);
    sizr = size(r,1);
    for jj = 1:sizr
        [r1,c1,v1] = find(TreeTab(:,2)==BrchD(r(jj),3));
         m = size(r1,1);
         if m == 0
        TreeTab(si,co) = BrchD(r(jj),3);
        TreeTab(si,co+BrchMax) = 1;
        TreeTab(si,co+2*BrchMax) = BrchD(r(jj),1);
        co = co+1;
        N = N+1;
         end
    end

    [r,c,v] = find(BrchD(:,3)==CurBus);
    sizr = size(r,1);
    for jj = 1:sizr
         [r1,c1,v1] = find(TreeTab(:,2)==BrchD(r(jj),2));
         m=size(r1);
             if m(1,1)==0
            TreeTab(si,co) = BrchD(r(jj),2);
            TreeTab(si,co+BrchMax) = 1;
            TreeTab(si,co+2*BrchMax) = BrchD(r(jj),1);
            co = co+1;
             N = N+1;
             end
    end

    % finding the nodes from the Trf
    [r,c,v] = find(TrfD(:,3)==CurBus);
    sizr = size(r,1);
    for jj = 1:sizr
         [r1,c1,v1] = find(TreeTab(:,2)==TrfD(r(jj),4));
         m=size(r1);
         if m(1,1)==0
        TreeTab(si,co) = TrfD(r(jj),4);
        TreeTab(si,co+BrchMax) = 2;
        TreeTab(si,co+2*BrchMax) = TrfD(r(jj),1);
        co = co+1;
         N = N+1;
         end
    end

    [r,c,v] = find(TrfD(:,4)==CurBus);
    sizr = size(r,1);
    for jj = 1:sizr
       [r1,c1,v1] = find(TreeTab(:,2)==TrfD(r(jj),3));
       m=size(r1);
             if m(1,1)==0  
            TreeTab(si,co) = TrfD(r(jj),3);
            TreeTab(si,co+BrchMax) = 2;
            TreeTab(si,co+2*BrchMax) = TrfD(r(jj),1);
            co = co+1;
             N = N+1;
             end
    end

    % finding the nodes from the Voltage Regulator
    if RegNo ~= 0 %&& RegD(1,1) ~=0
    [r,c,v] = find(RegD(:,3)==CurBus);
    sizr = size(r,1);
            for jj = 1:sizr

                [r1,c1,v1] = find(TreeTab(:,2)==RegD(r(jj),2));
                 m=size(r1);
                 if m(1,1)==0
                TreeTab(si,co) = RegD(r(jj),2);
                TreeTab(si,co+BrchMax) = 3;
                TreeTab(si,co+2*BrchMax) = RegD(r(jj),1);
                co = co+1;
                 N = N+1;
                 end
            end

            [r,c,v] = find(RegD(:,2)==CurBus);
            sizr = size(r,1);
            for jj = 1:sizr
                [r1,c1,v1] = find(TreeTab(:,2)==RegD(r(jj),3));
                 m=size(r1)
                     if m(1,1)==0  
                    TreeTab(si,co) = RegD(r(jj),3);
                    TreeTab(si,co+BrchMax) = 3;
                    TreeTab(si,co+2*BrchMax) = RegD(r(jj),1);
                    co = co+1;
                     N = N+1;
                     end
            end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % start of the counter logic for storing the buses to analyze virtually
    % like a pointer in memory

    if size(Temp,2) ~= 1 %&& Temp(1,1)==0
    [r1,c1,v1]= find(Temp);
    a = size(c1);
    b = a(1,2);
    if b~=0
    mam = c1(1,b);
    end
    else
        mam = 0;
    end

     for ii = 1:N
              Temp(1,mam+ii)= TreeTab(si,3+ii);
     end

    % Terminating the algorithm if the analysis of all nodes are done
    [r,c,v] = find(Temp);

        if size(c,2)==0 %&& CurBus ==0

        disp(' End of Tree forming algorithm');
            break;
        end

    % finding the non zero element in the 'Temp' matrix and ll - gives the 
    % column no of the last element in the matrix
    % c(1,ll) - means the last element of the Temp matrix

    ll = size(c,2);
    CurBus = Temp(1,c(1,ll));
    [r1,c1,v1] = find(Temp(1,:)==CurBus);

    for iii = 1:size(r1,2)                
               [r1,c1,v1] = find(Temp(1,:)== CurBus);            
                if size(c1,2)~=0
                    for zz = 1:size(c1,2)
                        Temp(1,c1(1,size(c1,2)))=0;
                    end
                end
    end

    NoTemp = size(Temp,2);
    NoTemp = NoTemp-1;
    si = si+1;

    end

    %TreeTab
    clk2 = clock;
    
    
    
%     save ('TreeTab.mat','TreeTab','BusChk','src','BrchMax','TermBus','clk2')
    save ('TreeTab.mat','TreeTab','src','BrchMax','TermBus')

    %TreeTab = [SiNo BusNo,BrchNo, RcvBus1, RcvBus2, ...... ,RcvBusN, Typ1, Typ2..
    %....,TypN, ElemId1, ElemId2, ElemId2,......ElemIdN] 

end
