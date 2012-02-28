%
% A function to generate a smooth beach profile on a uniform grid.  Takes
% the parameters interactively.  
%
function BeachGen
L0=input('Flat bottom run-up length at the left end L_0?');
L1=input('Sloping beach length L_1=?'); h0=input('h0?'); h1=input('h1?');
n1=input('Number of grid points on the sloping beach part n_1=?'); 
dx=L1/(n1-1); n0=round(L0/dx); L0=n0*dx; n=n0+n1; L=L0+L1;
display(['Run-up length re-defined to be L_1=',num2str(L1)]);
dh=h0-h1; x=[0:n-1]*dx; x0=L1/2;
for i=1:n0
    h(i)=h0;
end
for i=1:n1
    xbar=-1+((i-1)*dx/x0);
    h(n0+i)=h0-.5*dh*(1+f(xbar));
    hx(n0+i)=dh/x0*fx(xbar);
    hxx(n0+i)=-dh/(2*x0^2)*fxx(xbar); 
end
figure(1);clf;
plot(x,h,'-o');
figure(2); clf;
plot(x,hx,'-+',x,hxx,'-*');
%
% Save the bathymetry in an output file for OceanWave3D.
%
Header=['% Smooth beach with nx=',num2str(n),...
    ' L0=',num2str(L0),' Lbeach=',num2str(L1),' h_0=',num2str(h0), ...
    ' and h_1=',num2str(h1)];
fid = fopen(['beach_',num2str(n)],'wb'); % Open file
fprintf(fid,'%s\n',Header); % Write header
fprintf(fid,'%d\n',0);
for i=1:n
    fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e\n',h(i),hx(i),hxx(i),0,0); 
end
fclose('all');
%
%
function z=f(x)
%
% Normalized smooth beach function on the interval -1:1
%
if x>=1
    z=1;
elseif x<=-1
    z=-1;
else
    z=tanh(sin(pi/2*x)/(1-x^2));
end
function z=fx(x)
%
% Normalized smooth beach function on the interval -1:1
%
if x>=1
    z=0;
elseif x<=-1
    z=0;
else
    z=sech(sin(pi*x/2)/(1-x^2))^2 ...
        *(pi*(-1+x^2)*cos(pi*x/2)-4*x*sin(pi*x/2))/(4*(-1+x^2)^2);
end

function z=fxx(x)
%
% Normalized smooth beach function on the interval -1:1
%
if x>=1
    z=0;
elseif x<=-1
    z=0;
else
    z=sech(sin(pi*x/2)/(1-x^2))^2*( ...
        2*pi*x*cos(pi*x/2)/(-1+x^2)^2 ...
        -8*x^2*sin(pi*x/2)/(-1+x^2)^3 + 2*sin(pi*x/2)/(-1+x^2)^2 ...
        +pi^2*sin(pi*x/2)/(-4+4*x^2) ...
        -(pi*(-1+x^2)*cos(pi*x/2)-4*x*sin(pi*x/2))^2*f(x) ...
        /(2*(-1+x^2)^4));
end

