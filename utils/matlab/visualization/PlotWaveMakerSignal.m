function [t,y,u]=PlotWaveMakerSignal(fname)
%
% Read in the flux boundary condition and plot it up for visualization.
% 
% Note that the read statements here are more sensitive than the Fortran
% statements in the OW3D code, so no additional characters should appear
% after the header line. 
%
% Written by H.B. Bingham (hbb@mek.dtu.dk) Feb. 2016.
%
fid = fopen(fname,'r');
header=fgetl(fid);
tmp=fscanf(fid,'%f %d %d \n',[3,1]); dt=tmp(1); nt=tmp(2); ny=tmp(3);
y=zeros(ny,1); u=zeros(nt,ny);

y=fscanf(fid,'%f',ny);

for it=1:nt
    tmp=fscanf(fid,'%f',ny+1);
    t(it)=tmp(1);
    u(it,:)=tmp(2:ny+1);
end

fclose(fid);
%
% Plot the flux in time and space.
%
figure
title('Flux at the Western boundary');
if ny==1
    plot(t,u)
    xlabel('t'); ylabel('u');
else
    [Y,T]=meshgrid(y,t);
    mesh(T,Y,u);
    xlabel('t'); ylabel('y'); zlabel('u');
end
end
