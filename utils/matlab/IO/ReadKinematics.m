%
% This script reads the unformatted binary kinematics output file from the
% OceanWave3D code.  
%
% Modify the number of bits and the file name to correspond to your system.
%  
%clear all; %close all
%
% Set Nbits=32 or 64 and set idn= kinematics file # to read before calling
% this script
%
if exist('Nbits')==0 | exist('idn')==0
    error('Set "Nbits"=32 or 64 and "idn"= kinematics file number to read');
end
%
% Turn this flag on to compute time-derivatives and pressures. 
%
% Pressure='yes';
Pressure='no';
%
% Turn this flag on to plot up the data.  What is plotted is controlled
% after the reading block.  
%
% Plots='yes';
Plots='no';
%
% Open the file
if idn<10
    fname=['Kinematics0',num2str(idn),'.bin'];
else
    fname=['Kinematics',num2str(idn),'.bin'];
end
fid1 = fopen(fname); % File ID number for kinematics data

% Check for problems
if fid1 == -1
    error(['Could not open file ',fname,', returned fid1 = - 1']);
end

% Choose number of bits
if Nbits == 32
    int_nbit = 'int';
    njunkreak = 2;
elseif Nbits == 64
    int_nbit = 'int64';
    njunkread = 2;
else
    disp(['Illegal value for Nbits, Nbits = ' int2str(Nbits)]);
end
%
% Read the data from the file
% These read statements must correspond exactly to what appears in the
% Fortran subroutine: <top dir>/src/IO/StoreKinematicData.f90
%
junk = fread(fid1,1,int_nbit); % Read as a junk parameter
xbeg = fread(fid1,1,'int'); %
xend = fread(fid1,1,'int'); %
xstride = fread(fid1,1,'int'); %
ybeg = fread(fid1,1,'int'); %
yend = fread(fid1,1,'int'); %
ystride = fread(fid1,1,'int'); %
tbeg = fread(fid1,1,'int'); %
tend = fread(fid1,1,'int'); %
tstride = fread(fid1,1,'int'); %
dt = fread(fid1,1,'double'); % Time step size
nz = fread(fid1,1,'int'); %

junk = fread(fid1,2,int_nbit); % Junk read statements are necessary for eol markers

nx=floor((xend-xbeg)/xstride)+1; ny=floor((yend-ybeg)/ystride)+1; nt=floor((tend-tbeg)/tstride)+1;
%
% A scratch vector for reading the data
%
tmp=zeros(nx*ny*max(nz,5),1);
%
% The x-y grid, the depth and bottom gradients for this slice of data
%
tmp(1:5*nx*ny)=fread(fid1,5*nx*ny,'double');
junk = fread(fid1,2,int_nbit);
%
x=zeros(nx,ny); x(:)=tmp(1:5:5*nx*ny);
y=zeros(nx,ny); y(:)=tmp(2:5:5*nx*ny);
h=zeros(nx,ny); h(:)=tmp(3:5:5*nx*ny);
hx=zeros(nx,ny); hx(:)=tmp(4:5:5*nx*ny);
hy=zeros(nx,ny); hy(:)=tmp(5:5:5*nx*ny);
%
% The sigma coordinate
%
for i=1:nz
    sigma(i)=fread(fid1,1,'double');
end
junk = fread(fid1,2,int_nbit);
%
% Initialize arrays for the solution on this slice
%
eta=zeros(nt,nx,ny); etax=zeros(nt,nx,ny); etay=zeros(nt,nx,ny); 
phi=zeros(nt,nz,nx,ny); w=zeros(nt,nz,nx,ny); 
u=zeros(nt,nz,nx,ny); uz=zeros(nt,nz,nx,ny);
v=zeros(nt,nz,nx,ny); vz=zeros(nt,nz,nx,ny); wz=zeros(nt,nz,nx,ny);
t=[0:nt-1]*dt*tstride;   % The time axis
%
% Read in the solution variables eta, gradeta, phi, u, v, w, dudz, dvdz.  
%
for it=1:nt
    try
        tmp(1:nx*ny)=fread(fid1,nx*ny,'double');
        eta(it,:)=tmp(1:nx*ny);
        junk = fread(fid1,2,int_nbit);
    catch
        warning(['Read failed at time step ',num2str(it)]);
        break
    end
    %
    tmp(1:nx*ny)=fread(fid1,nx*ny,'double');
    etax(it,:)=tmp(1:nx*ny);
    junk = fread(fid1,2,int_nbit);
    %
    tmp(1:nx*ny)=fread(fid1,nx*ny,'double');
    etay(it,:)=tmp(1:nx*ny);
    junk = fread(fid1,2,int_nbit);
    %
    tmp=fread(fid1,nx*ny*nz,'double');
    phi(it,:)=tmp;
    junk = fread(fid1,2,int_nbit);
    %
    tmp=fread(fid1,nx*ny*nz,'double');
    u(it,:)=tmp;
    junk = fread(fid1,2,int_nbit);
    %
    tmp=fread(fid1,nx*ny*nz,'double');
    v(it,:)=tmp;
    junk = fread(fid1,2,int_nbit);
    %
    tmp=fread(fid1,nx*ny*nz,'double');
    w(it,:)=tmp;
    junk = fread(fid1,2,int_nbit);
    %
    [tmp,count]=fread(fid1,nx*ny*nz,'double');
    % Check for an incomplete run.
    if count==0, it=it-1; break; end
    wz(it,:)=tmp;
    junk = fread(fid1,2,int_nbit);
    %
    [tmp,count]=fread(fid1,nx*ny*nz,'double');
    % Check for an incomplete run.
    if count==0, it=it-1; break; end
    uz(it,:)=tmp;
    junk = fread(fid1,2,int_nbit);
    %
    [tmp,count]=fread(fid1,nx*ny*nz,'double');
    % Check for an incomplete run.
    if count==0, it=it-1; break; end
    vz(it,:)=tmp;
    junk = fread(fid1,2,int_nbit);
end
display(['Read ',num2str(it),' data points out of ',num2str(nt)])
%%
switch Pressure
    case 'yes'
        %
        % Compute the pressure and acceleration from the standard output
        % kinematics.  This is only done along the first slice in y for 3D
        % problems.
        %
        % Build the 4th-order even grid time differentiation matrix
        %
        alpha=2; r=2*alpha+1; c=BuildStencilEven(alpha,1);
        Dt=spdiags([ones(nt,1)*c(:,alpha+1)'],[-alpha:alpha],nt,nt);
        for j=1:alpha
            Dt(j,:)=0; Dt(j,1:r)=c(:,j)';
            Dt(nt-j+1,:)=0; Dt(nt-j+1,nt-r+1:nt)=c(:,r-j+1)';
        end
        Dt=Dt/dt;
        %
        % Compute time-derivatives of eta, phi, and u and eta_tt
        %
        % ip=input('x point index to work with?');
        etat=zeros(nt,nx); etatt=etat; phit=zeros(nt,nz,nx); p=phit; ut=p;
        for ip=1:nx
            etat(:,ip)=Dt*eta(:,ip); etatt(:,ip)=Dt*etat(:,ip);
            %
            for j=1:nz
                phit(:,j,ip)=Dt*phi(:,j,ip)-w(:,j,ip).*sigma(j).*etat(:,ip);
                p(:,j,ip)=-(phit(:,j,ip)+1/2*(u(:,j,ip).^2+w(:,j,ip).^2));
                ut(:,j,ip)=Dt*u(:,j,ip)-uz(:,j,ip).*sigma(j).*etat(:,ip);
            end
        end
end
%%
% A block for plotting the results
%
switch Plots
    case 'yes'
        %
        % ip=round(nx/2)+1; jp=1; % The horizontal grid point position to plot out
        % itp=[nt-33:2:nt]; % The time slice to focus on
        ip=1; jp=1; % The horizontal grid point position to plot out
        itp=[it-20:it]; % The time slice to focus on
        %
        % Plot the elevation as a function of time, then space at the last time
        % step.
        %
        ifig=1; figure(ifig);clf;
        plot(t,eta(:,ip,jp)); xlabel('t'); 
        ylabel(['eta(t) at x=',num2str(x(ip,jp))]);
        ifig=ifig+1; figure(ifig);clf;
        plot(x(:,1),eta(end,:,1),'-o',x(:,1),etax(end,:,1),'-+');
        xlabel('x'); title(['Free surface at t=',num2str(t(end))]);
        legend('eta(x)','eta_x(x)');
        %
        % This is the vertical coordinate for a linear problem: 
        %
        z=sigma*h(ip,jp)-h(ip,jp);
        %
        %
        ifig=ifig+1; figure(ifig);clf;
        plot(phi(itp,:,ip,jp),sigma,'-+');
        title(['Profiles of phi(sigma) at x=',num2str(x(ip,jp))] );
        ylabel('sigma');
        %
        ifig=ifig+1; figure(ifig);clf;
        plot(u(itp,:,ip,jp),sigma,'-+');
        title(['Profiles of u(sigma) at x=',num2str(x(ip,jp))] );
        ylabel('\sigma');
        %
        ifig=ifig+1; figure(ifig);clf;
        plot(w(itp,:,ip,jp),sigma,'-+');
        title(['Profiles of w(sigma) at x=',num2str(x(ip,jp))] );
        ylabel('sigma');
        %
        ifig=ifig+1; figure(ifig);clf;
        plot(uz(itp,:,ip,jp),sigma,'-+');
        title(['Profiles of du/dz(sigma) at x=',num2str(x(ip,jp))] );
        ylabel('sigma');
end
%