function wavegroupBredmoseJFM2009
close all
%clear all
clc

nx=1025;
H  = 1.45;
Lx = 50;
x0 = -200; % shifted to match a wave crest for this setup
h  = 4.25;
g  = 9.81;

SFsol.L       = Lx;
SFsol.k       = 2*pi/SFsol.L;
SFsol.h       = h;
SFsol.g       = g;
SFsol.uEorS   = 0;
SFsol.nsteps  = 8;
SFsol.maxiter = 64;
SFsol.n       = 32;
SFsol.EorS    = 0;
SFsol.H       = H;

[c,ubar,R,Q,cofA,cofB] = ExactSFSolution(SFsol.n,SFsol.H,SFsol.h,SFsol.L,SFsol.g,...
    SFsol.uEorS,SFsol.EorS,SFsol.nsteps,SFsol.maxiter);

x = linspace(-400,70,nx);
t = 0;

[U, W, E, Ex, P, Exx] = streamfunctionsol1D(SFsol.H,c,SFsol.k,SFsol.h,ubar,t,x,cofA,cofB);

modulation = sech(SFsol.k*(x'-x0)/4);
E=modulation.*E; P=modulation.*P;
plot(x,E,'k')


%
% Open a file for the initial condition to OceanWave3D
%
fid = fopen('OceanWave3D.init','w'); 
% 
% Write the initital condition to the file.  
%
fprintf(fid,'%s\n','% Initial condition for the Bredmose packet'); % Write header
fprintf(fid,'%f %f %d %d %f \n',[470,0,nx,1,0]);
for j=1:nx
    fprintf(fid, '%e %e \n',[E(j,1),P(j,1)]);
end


axis([-400 70 -4 2])
set(gcf,'position',[560         632        1010         294])
return

function [c,ubar,R,Q,A,B] = ExactSFSolution(nn,H,hh,LL,gg,uEorSS,EorSS,nsteps,maxiter)
%
% REFERENCE: Fourier method by Fenton & Rienecker (1982)
%
% Implemented by Allan P. Engsig-Karup, apek@imm.dtu.dk.
global n g k h Hi EorS uEorS
n = nn;
g = gg;
L = LL;
h = hh;
uEorS = uEorSS;
EorS  = EorSS;

% initialize some constants
k     = 2*pi/L;
T0    = 2*pi/sqrt(g*k*tanh(k*h));

% The first step in height and the solution from linear theory
Hi    = H/nsteps;
c     = L/T0;
ubar  = L/T0;
R     = g/(2*k)*tanh(k*h);
Q     = 0;
eta   = Hi/2*cos((1:n+1)*L/(2*n));
B     = zeros(1,n);
B(1)  = g*Hi/2*T0/L;

% The first solution step
inits = [ c;ubar;R;Q;eta(:);B(:)];
options = optimset('Display','iter','MaxIter',maxiter,'TolFun',1e-16,'Display','off');
f1    = fsolve(@solver,inits,options);
f2    = f1;%fsolve(@solver,f1,options);  % necessary??
for is = 2 : nsteps
    Hi = is*H/nsteps;
    inits = 2*f2-f1;
    f = fsolve(@solver,inits,options);
    f1 = f2;
    f2 = f;
end
% extract data
c    = f2(1);
ubar = f2(2);
R    = f2(3);
Q    = f2(4);
eta  = f2(5:5+n)';
B    = f2(5+n+1:5+2*n);
for i = 1 : n
    A(i) = 2/n*(1/2*(eta(1)+eta(n+1)*cos(k*i*L/2)));
    for j = 2 : n
        A(i) = A(i) + 2/n*eta(j)*cos((j-1)*k*i*L/(2*n));
    end
end

function f = solver(param)
global n g k h Hi EorS uEorS
c    = param(1);
ubar = param(2);
R    = param(3);
Q    = param(4);
eta  = param(5:5+n)';
B    = param(5+n+1:5+2*n);

% Compute estimates
i   = 1:n+1;
j   = 1:n;
psi = -ubar*eta + (((B'./(j*k))./cosh(j*k*h)))*(sinh(j'*k*(eta+h)).*cos(j'*(i-1)*pi/n)) ;
u   = -ubar + B'./cosh(j*k*h)*(cosh(j'*k*(eta+h)).*cos(j'*(i-1)*pi/n));
w   = B'./cosh(j*k*h)*(sinh(j'*k*(eta+h)).*sin(j'*(i-1)*pi/n));

% Setup equations
eq1 = eta(1) - eta(n+1) - Hi;
eq2 = eta(1) + eta(n+1) + 2*sum(eta(2:n));
if EorS==0
    eq3 = uEorS + ubar - c;
else
    if k*h<24
        eq3 = (uEorS + ubar - c)*h-Q;
    else
        eq3 = (uEorS + ubar - c)*h;
    end
end
bc1 = psi - Q;
bc2 = g*eta + 1/2*(u.^2+w.^2)-R;

% Return values
f = [ eq1(:) ; eq2(:) ; eq3(:); bc1(:); bc2(:) ];
return

function [U, W, E, Ex, P, Exx] = streamfunctionsol1D(H,c,k,h,ubar,t,X,cofA,cofB)
% Compute nonlinear stream function SURFACE solution in 1D
% SURFACE SOLUTION is given at z = eta.
%
% TO TEST SCRIPT: 
%   Make sure that derivatives of P matches U and W.
%
% By Allan P. Engsig-Karup, 02-02-2005.

if length(h)==1
    h = h*ones(size(X));
end

x = X(:)-c*t;

% Surface elevation
j   = 1 : length(cofA);
E   = cos(k*x*j)*cofA(:);
Ex  = -k*sin(k*x*j)*(cofA(:).*j');
Exx = -k^2*cos(k*x*j)*(cofA(:).*(j.^2)');

% Velocity components
% z  = E(:)*0-h(:); % check bottom
z  = E(:); % check bottom
U  = (c - ubar)   + ( cosh( k*(z(:)+h(:))*j )./cosh(k*h(:)*j ) .* cos(k*x*j) )*cofB(:);
W  = ( sinh( k*(z(:)+h(:))*j )./cosh(k*h(:)*j ) .* sin(k*x*j) )*cofB(:);

% Velocity Potential function
P  = (c - ubar)*x + ( cosh( k*(z(:)+h(:))*j )./cosh(k*h(:)*j ) .* sin(k*x*j) )*(cofB(:)./(j(:)*k));
% Px = U
return

