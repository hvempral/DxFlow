% The following code generates a report  
%
%               Programmer: Hemanth Kumar V, Michigan Technological Univ
%               Advisor: Dr Sumit Paudyal, MTU
%               Last Modified: 29th Jan 2015


clc
% clear all
global fullpathstr
load LFSresult.mat
disp(' Start of report writing')
D = [1 -1 0; 0 1 -1; -1 0 1];
clk = clock;

%% Print settings
inp_print = 00;
print_busV = 1;
print_busI = 0;
print_FeederI = 1;
print_trafoI = 1;


%%
% test = '\\homedir.mtu.edu\home\Desktop\DPFS\results';

    
file1 = [out filesep];
[pathstr,name,ext] = fileparts(filename);
filename1 = strcat(name,'.txt');
fullpathstr = strcat(file1,filename1);
fileID = fopen (fullpathstr,'w');
fprintf(fileID, '                        Load Flow Results  (Outfile Generator)                    \r\n');
fprintf(fileID, 'Report Generated on %s\r\n\r\n ',datestr(clk));
fprintf(fileID,'Load flow iteration: %d\n\r\n',iter);

fprintf(fileID,'-------------------------------------------------------------------------\r\n\r\n');

%fprintf(fileID, [repmat('\n %f\t %f\t %f\t\r', 1, size(t_a, 2)) '\n'], t_a')

if (inp_print==1) % print only when this flag is equal to 1
% input file generator
x = '%d\t %d\t\t %d\t %.3f\t\t %.3f\t %.3f\t \r\n\n'
fprintf(fileID, 'Bus Data:\r\n')
fprintf(fileID, 'BusNo\t BusName\t ZoneNo\t Voltage\t Vmin\t Vmax\t\r\n')
fprintf(fileID, '-----\t -------\t ------\t -------\t -----\t -----\t\r\n')
fprintf(fileID, x, (BusD)')

x = '%d\t %d\t\t %d\t %.3f\t %d\t %d\t \r\n\n'
fprintf(fileID, '\r\n Branch Data:\r\n')
fprintf(fileID, 'BrcNo\t FromNode\t ToNode\t length\t\t Config\t Status\t\r\n')
fprintf(fileID, '-----\t -------\t ------\t -------\t -----\t -----\t\r\n')
fprintf(fileID, x, (BrchD)')

x = '%d\t \r%d\t %d\t %.2f\t %.2f\t\t %.2f\t %.2f\t\t %.2f\t %.2f\t\t %d\t\t %d\t\t %d\t\t %d\t\t %d\t  \n\n'
fprintf(fileID, '\r\n Load Data:\r\n')
fprintf(fileID, 'Node\t Type\t Model\t Ph1(kW)\t Ph1(kVAR)\t Ph2(kW)\t  Ph2(kVAR)\t Ph3(kW)\t  Ph3(kVAR)\t Ph1-flg\t   Ph2-flg\t Ph3-flg\t LoadData\t  Status\t\r\n')
fprintf(fileID, '----\t ----\t -----\t -------\t ---------\t --------\t --------\t  --------\t --------\t  --------\t --------\t --------\t --------\t --------\t               \r\n')
fprintf(fileID, x, (LoadsD)')


 end
 
 
 
if (print_busV==1.0)
% fprintf(fileID,'\r\n\t\t\t ********************************************************************** \r\n');
fprintf(fileID, 'Voltage output: Node - wise\r\n');
% fprintf(fileID,'\r\n\t\t\t ********************************************************************** \r\n');
for ii = 1:length(BusD(:,1))
    vll(:,:,BusD(ii,1)) = D*v(:,:,BusD(ii,1));
    Vlnp(:,:,BusD(ii,1)) = W*vll(:,:,BusD(ii,1));
    Vllp(:,:,BusD(ii,1)) = D*Vlnp(:,:,BusD(ii,1));
    fprintf(fileID,'\r\nBus Node: %4.0f \r\n',BusD(ii,1));
    fprintf(fileID,'Voltage(V_ln)\t  Angle(deg) \t Voltage (V_LL)\t  Angle(deg) \r\n') ;
    for jj = 1:3
    fprintf(fileID,'%6.4f    @ %7.4f \t\t  %6.4f    @ %7.4f  \r\n', abs(v(jj,:,BusD(ii,1))),angle(v(jj,1,BusD(ii,1)))*180/pi,abs(vll(jj,:,BusD(ii,1))),angle(vll(jj,1,BusD(ii,1)))*180/pi);%,abs(W*v(jj,:,BusD(ii,1))),angle(W*v(jj,1,BusD(ii,1)))*180/pi);
%     fprintf(fileID,'%6.2f    @ %7.2f \t\t  %6.2f    @ %7.2f \t\t %6.2f    @ %7.2f \t\t  %6.2f    @ %7.2f\r\n', abs(v(jj,:,BusD(ii,1))),angle(v(jj,1,BusD(ii,1)))*180/pi, abs(vll(jj,:,BusD(ii,1))),angle(vll(jj,1,BusD(ii,1)))*180/pi, abs(Vlnp(jj,:,BusD(ii,1))),angle(Vlnp(jj,1,BusD(ii,1)))*180/pi, abs(Vllp(jj,:,BusD(ii,1))),angle(Vllp(jj,1,BusD(ii,1)))*180/pi);%,abs(W*v(jj,:,BusD(ii,1))),angle(W*v(jj,1,BusD(ii,1)))*180/pi);

    %fprintf(fileID,'%4.3f\t   @ %3.3f \t\t   %4.3f\t   @ %3.3f \r', abs(Vm_4(ii,1)),angle(Vm_4(ii,1))*180/pi,abs(VmLL_4(ii,1)),angle(VmLL_4(ii,1))*180/pi)
    end
    fprintf(fileID,'\n\n\n');
    fprintf(fileID,'-------------------------------------------------------------------------\r\n\r\n');
end
end


if (print_busI==1.0)
% fprintf(fileID,'\r\n\t\t\t ********************************************************************** \r\n');
fprintf(fileID, 'Current output: Node - wise:\r\n');
% fprintf(fileID,'\r\n\t\t\t ********************************************************************** \r\n');

for ii = 1:length(BusD(:,1))
    fprintf(fileID,'Current flowing out of Node: %4.0f \r\n',BusD(ii,1));
    fprintf(fileID,'Current(A)\t  angle(deg) \r\n');
    for jj = 1:3
    fprintf(fileID,'%4.5f      @ %10.5f \r\n', abs(Ibus(jj,:,BusD(ii,1))),angle(Ibus(jj,:,BusD(ii,1)))*180/pi);
    end
    fprintf(fileID,'-------------------------------------------------------------\r\n');
end
end



if (print_FeederI==1.0)
fprintf(fileID, '\r\n Feeder Currents: \r\n');
for ii = 1:length(BrchD(:,1))
    fprintf(fileID,'Current flowing from side of feeder: %4.0f [ Node %d - Node %d ] \r\n', BrchD(ii,1), BrchD(ii,2), BrchD(ii,3)) %,'(', ,',',BrchD(ii,3), ')');
    fprintf(fileID,'Current(A)\t  angle(deg) \r\n');
    for jj = 1:3
    fprintf(fileID,'%4.5f      @ %10.5f \r\n', abs(I(jj,:,BrchD(ii,1),1)),angle(I(jj,:,BrchD(ii,1),1))*180/pi);
    end
    fprintf(fileID,'-------------------------------------------------------------\r\n');
end
end


if (print_trafoI==1.0)
fprintf(fileID, '\r\n Transformer Currents: \r\n');
for ii = 1:length(TrfD(:,1))
    if (TrfD(ii,1)~= 0)
    fprintf(fileID,'Current flowing from side of Transformer: %4.0f [ Node %d - Node %d ] \r\n', TrfD(ii,1), TrfD(ii,3), TrfD(ii,4)) %,'(', ,',',BrchD(ii,3), ')');
    fprintf(fileID,'Current(A)\t  angle(deg) \r\n');
    for jj = 1:3
    fprintf(fileID,'%4.5f      @ %10.5f \r\n', abs(I(jj,:,TrfD(ii,1),2)),angle(I(jj,:,TrfD(ii,1),2))*180/pi);
    end
    end
    fprintf(fileID,'-------------------------------------------------------------\r\n');
end
end


fprintf(fileID,'===============================================================\r\n');


    if fileID~=1
          fclose(fileID);
    end
winopen(fullpathstr)
disp(' End of report writing')
