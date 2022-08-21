%
% Validation of the kinematics for a stream function wave.  This check is
% for the test case set up in OceanWave3D.inp.StreamFuncWave and compares
% the kinematics at peak u and w points with those computed directly from
% the function 'StreamFunctionCoefficients.m', which must be in the path. 
%
% The OW3D utility 'ReadKinematics.m' must also be in the path. 
%
% Written by H.B. Bingham, DTU, August 21, 2022
%
%%
%
clear all;
%
Nbits=32;  idn=1; ifig=0;
%
% Read the kinematics file from the OW3D run.
%
ReadKinematics
%
% Compare the calculations with the OceanWave3D.inp.StreamFuncWave case 
% for a Stream function wave H/L=0.08, kh=2pi, U_Stokes=0.
%
g=9.82; H=.08; L=1; k=2*pi/L; h0=1; kh=k*h0; nSF=24; nsteps=6; uEorS=0;
EorS='Stokes'; TorL='L';
%
% Compute all the Stream Function Theory parameters for this case. 
%
[cSF q R k etaN Bsf(:,1)]=StreamFunctionCoefficients(nSF,H,...
    h0,L,TorL,uEorS,EorS,nsteps,g); 
%
% Get the Fourier coefficients of the elevation so it can be evaluated at
% arbitrary x-positions.
%
js=[1:nSF-1]; dx=L/(2*nSF);
for i=1:nSF
    Asf(i,1)=2/nSF*(.5*(etaN(1)+etaN(nSF+1)*cos(k*i*L/2)) ...
        + sum(etaN(2:nSF).*cos(js*k*i*dx)));
end
% Multiply Nyquist by 1/2 so that straight sums can be
% used later.
%
Asf(nSF,1)=1/2*Asf(nSF,1);
%
% Form the full j-vector for computing the kinematics. 
%
js=[1:nSF]';
%
% The wave period and frequency.
%
T=L/cSF; om=2*pi/T;
%%
% These are the resolutions from the OceanWave3D.inp file. Pick out one 
% wavelength at the center of the domain at the last time step. 
%
ppw=16; ip0=(nx-1)/2-ppw/2+1; ip=[ip0:ip0+ppw]; 
it=nt;
time=t(it);
%
% Compute the surface values and compare.
%
for i=1:length(ip)
    eta2(i,1)=sum(Asf.*cos(js*k*(x(ip(i))-cSF*time)));
    phiS2(i,1)=((cSF-Bsf(1,1))*x(ip(i))...
        +sum(Bsf(2:nSF+1,1)./(js*k).*cosh(k*js.*(eta2(i,1)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x(ip(i))-cSF*time))));
    uS2(i,1)=((cSF-Bsf(1,1))...
        +sum(Bsf(2:nSF+1,1).*cosh(k*js.*(eta2(i,1)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x(ip(i))-cSF*time))));
    wS2(i,1)=sum(Bsf(2:nSF+1,1).*sinh(k*js.*(eta2(i,1)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x(ip(i))-cSF*time)));
end

ifig=ifig+1; figure(ifig);
plot(x(ip),eta(it,ip)/H,'-+',x(ip),eta2(:,1)/H,'o'); grid on;
xlabel('x/L'); ylabel('\eta/H')
legend('OW3D','SF')
title(['t=',num2str(t(it))]);

ifig=ifig+1; figure(ifig);
plot(x(ip),squeeze(phi(it,nz,ip)),'-+',x(ip),phiS2(:,1),'o'); grid on;
xlabel('x/L'); ylabel('\phi_S')
legend('OW3D','SF')
title(['t=',num2str(t(it))]);

ifig=ifig+1; figure(ifig);
plot(x(ip),squeeze(u(it,nz,ip)),'-+',x(ip),uS2(:,1),'o'); grid on;
xlabel('x/L'); ylabel('u_S')
legend('OW3D','SF')
title(['t=',num2str(t(it))]);

ifig=ifig+1; figure(ifig);
plot(x(ip),squeeze(w(it,nz,ip)),'-+',x(ip),wS2(:,1),'o'); grid on;
xlabel('x/L'); ylabel('w_S')
legend('OW3D','SF')
title(['t=',num2str(t(it))]);

%%
% Compare the kinematics under the wave crest.
%
ipc=ip0+ppw/2;
z=sigma*(h(ipc,1)+eta(it,ipc))-h(ipc,1); x0=x(ipc,1); eta0=eta2(ppw/2+1,1);
for i=2:length(z)
    phiSF(i-1,1)=((cSF-Bsf(1,1))*x0...
        +sum(Bsf(2:nSF+1,1)./(js*k).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time))));
    uSF(i-1,1)=((cSF-Bsf(1,1))...
        +sum(Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time))));
    uxSF(i-1,1)=(-sum(k*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time))));
    uzSF(i-1,1)=sum(k*js.*Bsf(2:nSF+1,1).*sinh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time)));
    wSF(i-1,1)=sum(Bsf(2:nSF+1,1).*sinh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time)));
    wxSF(i-1,1)=sum(k*js.*Bsf(2:nSF+1,1).*sinh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time)));
    wzSF(i-1,1)=sum(k*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time)));
    pSF(i-1,1)=-(-sum(cSF*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time))) + .5*(uSF(i-1,1).^2 ...
        + wSF(i-1,1).^2));
    % Lagrangian particle accelerations in the horizontal and vertical: 
    utSF(i-1,1)=sum(om*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time))) + uSF(i-1,1)*uxSF(i-1,1) ...
        + wSF(i-1,1)*uzSF(i-1,1);
    wtSF(i-1,1)=-sum(om*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time))) + uSF(i-1,1)*wxSF(i-1,1) ...
        + wSF(i-1,1)*wzSF(i-1,1);
end
%
% Potential and velocities
%
ifig=ifig+1; figure(ifig);
plot(phi(it,2:nz,ipc),z(2:nz),'-+',phiSF(:,1),z(2:nz),'o',...
    u(it,2:nz,ipc),z(2:nz),'-+',uSF(:,1),z(2:nz),'o',...
    w(it,2:nz,ipc),z(2:nz),'-+',wSF(:,1),z(2:nz),'o')
legend('\phi','\phi_e','u','u_e','w','w_e','Location','SouthEast');
ylabel('z/L');
title('Crest values'); grid on
%
% Gradient of the velocities
%
ifig=ifig+1; figure(ifig);
plot(uz(it,2:nz,ipc),z(2:nz),'-+',uzSF(:,1),z(2:nz),'o',...
    ux(it,2:nz,ipc),z(2:nz),'-+',uxSF(:,1),z(2:nz),'o',...
    wz(it,2:nz,ipc),z(2:nz),'-+',wzSF(:,1),z(2:nz),'o',...
    wx(it,2:nz,ipc),z(2:nz),'-+',wxSF(:,1),z(2:nz),'o');
legend('u_z','(u_z)_e','u_x','(u_x)_e','w_z','(w_z)_e','w_x','(w_x)_e','Location','SouthEast');
ylabel('z/L');
title('Crest values'); grid on
%
% Errors in the Laplacian of phi and cross derivatives.
%
ifig=ifig+1; figure(ifig);
plot(ux(it,2:nz,ipc)+wz(it,2:nz,ipc),z(2:nz),'-+',...
    uz(it,2:nz,ipc)-wx(it,2:nz,ipc),z(2:nz),'-+');
legend('Errors in u_x+w_z','Errors in u_z-w_x','Location','SouthEast');
ylabel('z/L');
title('Crest values'); grid on
%
% Pressure and Lagrangian particle acceleration
%
utL=ut+u.*ux+w.*uz; 
wtL=wt+u.*wx+w.*wz; 
%
ifig=ifig+1; figure(ifig);
plot(p(it,2:nz,ipc),z(2:nz),'-+',pSF(:,1),z(2:nz),'o',...
    utL(it,2:nz,ipc),z(2:nz),'-+',utSF(:,1),z(2:nz),'o',...
    wtL(it,2:nz,ipc),z(2:nz),'-+',wtSF(:,1),z(2:nz),'o')
legend('p','p_e','u_t','(u_t)_e','w_t','(w_t)_e','Location','SouthWest');
ylabel('z/L');
title('Crest values of p and Lagrangian accelerations'); grid on

%%
% Compare the kinematics near the peak vertical velocity. This is
% hard-coded to the 16 points per wavelength resolution.
%
ipc=ip0+ppw/2+4;
z=sigma*(h(ipc,1)+eta(it,ipc))-h(ipc,1); x0=x(ipc,1); eta0=eta2(ppw/2+5,1);
for i=2:length(z)
    phiSF(i-1,1)=((cSF-Bsf(1,1))*x0...
        +sum(Bsf(2:nSF+1,1)./(js*k).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time))));
    uSF(i-1,1)=((cSF-Bsf(1,1))...
        +sum(Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time))));
    uxSF(i-1,1)=(-sum(k*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time))));
    uzSF(i-1,1)=sum(k*js.*Bsf(2:nSF+1,1).*sinh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time)));
    wSF(i-1,1)=sum(Bsf(2:nSF+1,1).*sinh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time)));
    wxSF(i-1,1)=sum(k*js.*Bsf(2:nSF+1,1).*sinh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time)));
    wzSF(i-1,1)=sum(k*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time)));
    pSF(i-1,1)=-(-sum(cSF*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time))) + .5*(uSF(i-1,1).^2 ...
        + wSF(i-1,1).^2));
    % Lagrangian particle accelerations in the horizontal and vertical: 
    utSF(i-1,1)=sum(om*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*sin(k*js*(x0-cSF*time))) + uSF(i-1,1)*uxSF(i-1,1) ...
        + wSF(i-1,1)*uzSF(i-1,1);
    wtSF(i-1,1)=-sum(om*js.*Bsf(2:nSF+1,1).*cosh(k*js.*(z(i)+h0))...
        ./cosh(js*k*h0) .*cos(k*js*(x0-cSF*time))) + uSF(i-1,1)*wxSF(i-1,1) ...
        + wSF(i-1,1)*wzSF(i-1,1);
end

ifig=ifig+1; figure(ifig);
plot(phi(it,2:nz,ipc),z(2:nz),'-+',phiSF(:,1),z(2:nz),'o',...
    u(it,2:nz,ipc),z(2:nz),'-+',uSF(:,1),z(2:nz),'o',...
    w(it,2:nz,ipc),z(2:nz),'-+',wSF(:,1),z(2:nz),'o')
legend('\phi','\phi_e','u','u_e','w','w_e','Location','SouthEast');
ylabel('z/L');
title('Near the peak vertical vel.'); grid on

ifig=ifig+1; figure(ifig);
plot(uz(it,2:nz,ipc),z(2:nz),'-+',uzSF(:,1),z(2:nz),'o',...
    ux(it,2:nz,ipc),z(2:nz),'-+',uxSF(:,1),z(2:nz),'o',...
    wz(it,2:nz,ipc),z(2:nz),'-+',wzSF(:,1),z(2:nz),'o',...
    wx(it,2:nz,ipc),z(2:nz),'-+',wxSF(:,1),z(2:nz),'o');
legend('u_z','(u_z)_e','u_x','(u_x)_e','w_z','(w_z)_e','w_x','(w_x)_e','Location','SouthEast');

ylabel('z/L');
title('Near the peak vertical vel.'); grid on

ifig=ifig+1; figure(ifig);
plot(ux(it,2:nz,ipc)+wz(it,2:nz,ipc),z(2:nz),'-+',...
    uz(it,2:nz,ipc)-wx(it,2:nz,ipc),z(2:nz),'-+');
legend('Errors in u_x+w_z','Errors in u_z-w_x','Location','SouthEast');
ylabel('z/L');
title('Near the peak vertical vel.'); grid on
%
% Pressure and Lagrangian acceleration
%
utL=ut+u.*ux+w.*uz; 
wtL=wt+u.*wx+w.*wz; 
%
ifig=ifig+1; figure(ifig);
plot(p(it,2:nz,ipc),z(2:nz),'-+',pSF(:,1),z(2:nz),'o',...
    utL(it,2:nz,ipc),z(2:nz),'-+',utSF(:,1),z(2:nz),'o',...
    wtL(it,2:nz,ipc),z(2:nz),'-+',wtSF(:,1),z(2:nz),'o')
legend('p','p_e','u_t','(u_t)_e','w_t','(w_t)_e','Location','SouthEast');
ylabel('z/L');
title('p and Lagrangian acceleration near peak w'); grid on
