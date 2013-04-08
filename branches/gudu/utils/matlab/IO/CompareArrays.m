close all
clear all
clc

% compare computed and exact solution (IC=8)

IOmethod = 1; %0:binary ; 1:classical unformatted ; 2:unformatted ftn95

% CHOOSE BETWEEN byte storage format below (uncomment the one needed)
% byteorder = 'ieee-be'; % IEEE Big-Endian format 
byteorder = 'ieee-le'; % IEEE Little-Endian format

filename = '/Users/apek/OW3Dtest/R.bin';
R = LoadRealArray(filename,IOmethod, byteorder);

filename = '/Users/apek/OW3Dtest/P.bin';
P = LoadRealArray(filename,IOmethod, byteorder);
return

filename='../W.bin';
W = LoadRealArray(filename,IOmethod, byteorder);

filename='../WA.bin';
WA = LoadRealArray(filename,IOmethod, byteorder);

filename='../X.bin';
X = LoadRealArray(filename,IOmethod, byteorder);

filename='../Y.bin';
Y = LoadRealArray(filename,IOmethod, byteorder);

% filename='../Z.bin';
% Z = LoadRealArray(filename,IOmethod, byteorder);

subplot(1,3,1)
mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),WA(2:end-1,2:end-1))
title('Exact W')

subplot(1,3,2)
mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),W(2:end-1,2:end-1))
title('Computed W')

subplot(1,3,3)
mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),W(2:end-1,2:end-1)-WA(2:end-1,2:end-1))
title('Error')

set(gcf,'position',[254         475        1518         546])