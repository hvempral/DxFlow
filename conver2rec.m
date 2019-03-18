% Matlab code for polar to rectangular function, as a part of 
% Three Phase load flow Program 
%               Programmer: Hemanth Kumar V, Michigan Technological Univ
%               Advisor: Dr Sumit Paudyal, MTU
%               Last Modified: 29th Jan 2015
function res = conver2rec(mag, ang)
    [x,y]=pol2cart(ang*pi/180,mag);
    res = x+i*y;
end