function OW3D_LinearDispersionErrorsWith_kh
%%
%
% This utility function can be used to produce linear dispersion error
% plots as a function of kh for one value of nx points per wavelenth, one
% value of stencil 1/2-width alpha, and several values of vertical grid
% points nz. The vertical grid can be uniform or cosine-spaced towards the
% free-surface. 
%
% See Engsig-Karup, Bingham and Lindberg, JCP (2009) Sec. 4.3 for a
% detailed discussion of the theory. 
%
% Written by H. B. Bingham, DTU Mechanical Engineering. 
%
%%
%
% Set nx, alpha, and max. kh and vertical grid spacing. 
%
nx=8;             % Number of grid points per wavelength. 
alpha=3;          % Stencil 1/2-width, r=2*alpha+1. 1 <= alpha <= 4. 
p=2*alpha; r=p+1;
kh1=30;           % kh_max
Spacing='Cosine'; % Vertical spacing can be 'Uniform' or 'Cosine'. 
%
%%
% Solve for the dispersion errors and plot.
%
nkh=50; dkh=kh1/(nkh-1); kkh=[.01 [1:nkh-1]*dkh];
k=2*pi;
%
switch Spacing
    case 'Uniform'     % Uniform grid
        nnz=4; nzvec=[max(round(nx/2),2*alpha+1) nx 2*nx 4*nx];
        for iz=1:length(nzvec)
            nz0=nzvec(iz); nz=nz0+1;
            for ikh=1:length(kkh)
                kh=kkh(ikh); h=kh/k;
                dz=h/(nz0-1); z=-dz*[nz0-1:-1:0]; z=[z(1)-(z(2)-z(1)) z];
                
                w=SolveForWk2PiGhost(nx,nz,alpha,z);
                disp(ikh,iz)=w(nz)/k;
            end
            dispE=tanh(kkh);
        end
    case 'Cosine'     % Clustered grid
        nnz=4;
        if alpha==1
            nzvec=[25 50 75 150];
        elseif alpha==2
            nzvec=[9 12 17 max(35,2*nx+1)];
        elseif alpha==3
            nzvec=[7 9 12 17];
        else
            nzvec=[p p+2 p+4 2*p];
        end
        for iz=1:length(nzvec)
            nz0=nzvec(iz); nz=nz0+1;
            for ikh=1:length(kkh)
                kh=kkh(ikh); h=kh/k;
                th0=0;
                th1=pi/2; th=th1-th0; dth=th/(nz0-1); l0=sin(th0); l1=sin(th1);
                l=l1-l0; zv=h*(-1+(-l0+sin(th0+dth*[0:nz0-1]))/l);
                zv=[zv(1)-(zv(2)-zv(1)) zv];
                
                w=SolveForWk2PiGhost(nx,nz,alpha,zv);
                disp(ikh,iz)=w(nz)/k;
            end
            dispE=tanh(kkh);
        end
end
%
% Plot the errors in dispersion.
%
ifig=0;
ifig=ifig+1;figure(ifig); clf;
plot(kkh,dispE,'-k',kkh,disp(:,1),kkh,disp(:,2),...
    kkh,disp(:,3),kkh,disp(:,3))
legend('Exact',['N_z=',num2str(nzvec(1))],['N_z=',num2str(nzvec(2))],...
    ['N_z=',num2str(nzvec(3))],['N_z=',num2str(nzvec(4))],'Location','SouthWest');
xlabel('kh'); ylabel('Dispersion');grid on; axis([0 kkh(end) 0 1.1]);
switch Spacing
    case 'Uniform'
        title(['p=',num2str(2*alpha),' uniform grid, N_x=',num2str(nx)]);
    case 'Cosine'
        title(['p=',num2str(2*alpha),' clustered grid, N_x=',num2str(nx)]);
end
%
ifig=ifig+1; figure(ifig); clf; plot(kkh,disp(:,1)-dispE',':',...
    kkh,disp(:,2)-dispE','-.',...
    kkh,disp(:,3)-dispE','--',kkh,disp(:,4)-dispE','-')
axis([0 kkh(end) -.02 .02]);
xlabel('kh'); ylabel('Relative dispersion error');
legend(['N_z=',num2str(nzvec(1))],['N_z=',num2str(nzvec(2))],...
    ['N_z=',num2str(nzvec(3))],['N_z=',num2str(nzvec(4))],'Location','NorthWest'); grid on;

switch Spacing
    case 'Uniform'
        title(['r=',num2str(r),', uniform vertical grid, N_x=',num2str(nx)]);
    case 'Cosine'
        title(['r=',num2str(r),', clustered vertical grid, N_x=',num2str(nx)]);
end
%%
% Included functions
%

function [w phi]=SolveForWk2PiGhost(nx,nz,alpha,z)
%%
% This function solves the discrete linear dispersion operator on a uniform
% horizontal grid using nx grid points per wavelength and a possibly
% non-uniform vertical grid spacing given by z(1:nz). The stencil-size 
% of the finite-difference operators is r=2*alpha+1 and they are centered
% in the horizontal but off-centered towards the ends in the vertical. One
% ghost layer is assumed beneath the horizontal bottom to impose both the
% homogeneous Neumann condition and the Laplace equation there. 
%
% Only implemented up to alpha=4 at this point. 
%
% Written by H.B. Bingham, DTU Mechanical Engineering. 
%
rank=2*alpha+1;
%
% The horizontal second-derivative operator. 
%
if(alpha==1)
  coeffx=nx^2*(-2+2*cos(2*pi/nx));
elseif(alpha==2)
  coeffx=nx^2*(-5/2+8/3*cos(2*pi/nx)-1/6*cos(4*pi/nx));
elseif(alpha==3)
  coeffx=nx^2*(-49/18+3*cos(2*pi/nx)-3/10*cos(4*pi/nx)+1/45*cos(6*pi/nx));
elseif(alpha==4)
  coeffx=nx^2*(-205/72+16/5*cos(2*pi/nx)-2/5*cos(4*pi/nx)+16/315*cos(6*pi/nx)...
	       -1/280*cos(8*pi/nx));
else
  error('This solution is only implemented up to alpha=4')
end
%
% The vertical first- and second-derivative operators.
%
dfdz=BuildStencilVariable(nz,alpha,z,1);
dfdzz=BuildStencilVariable(nz,alpha,z,2);
%
% Build the matrix for this line of vertical grid points.
%
A=zeros(nz);
A(1,1:rank)=dfdz(1:rank,2); % Bottom condition at grid point 2.
for im=2:alpha              % Laplace equation
  A(im,1:rank)=dfdzz(1:rank,im);
end
for im=alpha+1:nz-alpha
  A(im,im-alpha:im+alpha)=dfdzz(1:rank,im);
end
for im=nz-alpha+1:nz-1
  A(im,nz-rank+1:nz)=dfdzz(1:rank,im);
end
A(nz,nz)=1;                 % Dirichlet condition at z=0
%
% Add in the second x-derivative operator.
%
for im=2:nz-1
  A(im,im)=A(im,im)+coeffx;
end
%
% Solver for phi and w.
%
rhs=zeros(nz,1); rhs(nz,1)=1;
phi=A\rhs;
w=DfDxVariable(nz,alpha,phi',dfdz);

function fx=BuildStencilVariable(N,alpha,x,der)
%%
% A function to compute finite-difference coefficients for the der'th 
% derivative of a function on a variably spaced grid x(1:N) using 
% 2*alpha+1 neighboring points.  N  sets of 
% coefficients are returned where set 1,2,...,alpha are one-sided schemes 
% at the left end, alpha+1,...,N-alpha are centered schemes, and 
% N-alpha+1,...,N  are one-sided schemes at the right end.  
% 
% Written by H.B. Bingham, DTU Mechanical Engineering. 

r=2*alpha+1;
%
% One-sided schemes for the left end-points.  
%
for ip=1:alpha
    c=FornbergInterpolate(x(ip),x(1:r),r,der); 
    fx(1:r,ip)=c(der+1,1:r)';
end
%
% The centered schemes
%
for ip=alpha+1:N-alpha
    c=FornbergInterpolate(x(ip),x(ip-alpha:ip+alpha),r,der);
    fx(1:r,ip)=c(der+1,1:r)';
end
%
% One-sided schemes for the right end-points.  
%
for ip=N-alpha+1:N
    c=FornbergInterpolate(x(ip),x(N-r+1:N),r,der); 
    fx(1:r,ip)=c(der+1,1:r)';
end


function c=FornbergInterpolate(x0,x,n,m)
%%
% A matlab implementation of the algorithm of 
%
%     B. Fornberg, 1998. "Calculation of weights in finite difference
%     formulas".  SIAM Rev. Vol. 40, No. 3, pp. 685-691.  
%
% for calculating arbitrary-order FD schemes on arbitrary 1-D grids for all
% applicable derivatives.  The interpolations are found at the point x0,
% using the function positions x(1:n), for derivatives 0:m, m<=n-1;   
%
% 
% Written by H.B. Bingham
%            Mech. Eng.
%            Technical U. of Denmark
%
% This version:  2009
%
%%
if m>n-1, error('Too many derivatives for this stencil, try again.'); end
%
c1=1; c4=x(1)-x0;
c=zeros(n,m+1);
c(1,1)=1;
for i=2:n
    mn=min(i,m+1); c2=1; c5=c4; c4=x(i)-x0;
    for j=1:i-1
        c3=x(i)-x(j);
        c2=c2*c3;
        if j==i-1
            for k=mn:-1:2
                c(k,i)=c1*((k-1)*c(k-1,i-1)-c5*c(k,i-1))/c2;
            end
            c(1,i)=-c1*c5*c(1,i-1)/c2;
        end
        for k=mn:-1:2
            c(k,j)=(c4*c(k,j)-(k-1)*c(k-1,j))/c3;
        end
        c(1,j)=c4*c(1,j)/c3;
    end
    c1=c2;
end
return

function fprime=DfDxVariable(N,alpha,f,fx)
%%
% This function takes the first derivative of the vector 
% f(1:N) using the N 2*alpha+1 point stencils in fx.  
% The stencils are one-sided at the ends and centered in 
% the middle. 
%
% Written by H. B. Bingham, DTU Mechanical Engineering. 
%
%
rank=2*alpha+1;

% The left end points
for ip=1:alpha
    fprime(ip)=f(1:rank)*fx(1:rank,ip);
end
% The mid-section
for ip=alpha+1:N-alpha
    fprime(ip)=f(ip-alpha:ip+alpha)*fx(1:rank,ip);
end  
%The right end points
for ip=N-alpha+1:N
    fprime(ip)=f(N-rank+1:N)*fx(1:rank,ip);
end
