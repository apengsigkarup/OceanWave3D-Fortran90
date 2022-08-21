function [c,q,R,k,eta,B,exitflag]=StreamFunctionCoefficients(n,H,h,...
                                    TorL,intype,uEorS,EorS,nsteps,g)
%
%  This function computes the stream function theory solution for a wave of
%  period T or wavelength L, height H, and mean Stokes drift/Eulerian 
%  velocity uEorS in water depth h using n modes and nsteps 
%  in nonlinearity.
%
%  intype='T' treats TorL as an input wave period T.
%  intype='L' treats TorL as an input wavelength L.
%
%  EorS='Euler' treats uEorS as the Eulerian mean velocity.
%  EorS='Stokes' treats uEorS as the Stokes drift velocity.
%
%  The returned solution is the celerity c, the flux q, the Bernoulli 
%  constant R, the wavenumber k, eta(1:n+1) the elevation at n+1 evenly 
%  spaced points from crest to trough, and the coefficients B(1:n+1)=B_j, 
%  j=0:n where
%
%  psi(x,z) = -B_0 z
%             + sum_{j=1}^n B_j/(jk) sinh(jk(z+h))/cosh(jkh) cos(jkx)
%
%  is the stream function.  The corresponding velocity components are
%  obtained by differentiating this expression:  
%
%      (u,w)=(dpsi/dz, -dpsi/dx).  
%  
%  The implementation follows the theory described in:  J.D. Fenton.   
%  The numerical solution of steady water wave problems. 
%  Comput. Geosci. 14:3, 357-68.  1988.
%
%  Written by Harry B. Bingham. Modified by Raphael Mounet to support both
%  wave length and period as input. 
%
%  This version:  20.04.2020

switch intype
    case 'T'
        T = TorL; omega = 2*pi/T; 
        syms k; k0 = eval(vpasolve(omega^2-g*k*tanh(k*h)==0,k,0.05));
        clear k;
        L0 = 2*pi/k0; c0 = L0/T;
        % 
        % The variables are ordered: x=[c,q,R,(eta_j,j=0:n),(B_j,j=0:n),k]
        % Provide the initial guess based on linear theory:
        %
        H0=H; H=H0/nsteps;
        x0=[c0 0 g/(2*k0)*tanh(k0*h) H/2*cos((0:n).*pi/n) ...
            c0 g*H/2*T/L0 zeros(1,n-1) k0];
        
    case 'L'
        L = TorL; k = 2*pi/L; 
        T0 = 2*pi/sqrt(g*k*tanh(k*h)); c0 = L/T0;
        % 
        % The variables are ordered: x=[c,q,R,(eta_j,j=0:n),(B_j,j=0:n)]
        % Provide the initial guess based on linear theory:
        %
        H0=H; H=H0/nsteps;
        x0=[c0 0 g/(2*k)*tanh(k*h) H/2*cos([0:n]*pi/n) ...
            c0 g*H/2*T0/L zeros(1,n-1)];
        
    otherwise
        error('intype must be either "T" or "L".');
end
%      
% Find the solution using a "paramatrized" function to work out the
% right hand sides of the stream function theory solution:
%
options=optimset('MaxFunEvals',5000,'MaxIter',2000);
[x,fval,exitflag]=fsolve(@(x) ...
    fStreamFunc(x,n,H,h,TorL,intype,uEorS,EorS,g),x0,options);
%
% Iterate up to the full height in nsteps, with the initial guess
% extrapolated from the previous solution.
%
for is=2:nsteps
    H=is*H0/nsteps; x0=2*x-x0;
    [x,fval,exitflag]=fsolve(@(x) ...
        fStreamFunc(x,n,H,h,TorL,intype,uEorS,EorS,g),x0,options);
end
%
% Extract the solution from the final solution vector:
%
c=x(1); q=x(2); R=x(3);
eta=x(4:4+n); B=x(n+5:2*n+5);
if length(x)==2*n+6
     k=x(end);
end

    function F=fStreamFunc(x,n,H,h,TorL,intype,uEorS,EorS,g)
        %
        % A function to evaluate the stream function theory right hand sides
        %
        % The variables are: 
        %    x = [c,q,R,(eta_j,j=0:n),(B_j,j=0:n),k] if k is unknown (case 
        %       with intype = 'T'), 
        % or x = [c,q,R,(eta_j,j=0:n),(B_j,j=0:n)] if k is already
        %       specified (case with intype = 'L')
        %
        % The equations are:
        %                    eta_0-eta_n-H=0
        %                    eta_0+eta_n+sum_{j=1}^{n-1} eta_j = 0
        %      for EorS='Euler'
        %                    -c + uEorS + B_0 = 0
        %      for EorS='Stokes'
        %                    (-c + uEorS + B_0)h = q
        %
        %  Equations 4:4+n are the kinematic conditions at each x: 
        %                    psi(eta) = q
        %
        %  Equations 4+n+1:2n+5 are the dynamic conditions at each x:
        %                    g eta + 1/2(u^2+w^2) = R
        %
        %  Equation 2n+6 is involved only if k is unknown and it 
        %  corresponds to the definition of the phase velocity c:
        %                    k-2*pi/(c*T) = 0
        %
        
        js=1:n;
        
        F(1)=x(4)-x(4+n)-H;
        F(2)=x(4)+x(4+n)+2*sum(x(5:4+n-1));
        
        switch EorS
            case 'Euler'
                F(3)=-x(1)+uEorS+x(n+5);
            case 'Stokes'
                F(3)=(-x(1)+uEorS+x(n+5))*h-x(2);
            otherwise
                error('EorS must be either "Euler" or "Stokes".');
        end
        
        switch intype
            case 'T'
                T = TorL;
                % Loop over the n-points from crest (x=0) to trough (x=L/2)
                for i=0:n
                    % The stream function and velocity at this point x:
                    psi=-x(n+5)*x(4+i)...
                        +sum(x(n+6:2*n+5)./(js*x(end)).*sinh(x(end)*js.*...
                        (x(4+i)+h))./cosh(js*x(end)*h).*cos(js*i*pi/n));
                    u=-x(n+5) ...
                        +sum(x(n+6:2*n+5).*cosh(x(end)*js.*(x(4+i)+h))...
                        ./cosh(js*x(end)*h).*cos(js*i*pi/n));
                    w= sum(x(n+6:2*n+5).*sinh(x(end)*js.*(x(4+i)+h))...
                        ./cosh(js*x(end)*h).*sin(js*i*pi/n));
                    
                    % The kinematic condition:  psi(eta)=q
                    F(4+i)=psi-x(2);
                    
                    % The dynamic condition:  g eta + 1/2(u^2+w^2)=R
                    F(4+n+1+i)=g*x(4+i) + 1/2*(u^2 + w^2) - x(3);
                end
                F(2*n+6)=x(end)-2*pi/(x(1)*T);
                
            case 'L'
                k = 2*pi/TorL;
                % Loop over the n-points from crest (x=0) to trough (x=L/2)
                for i=0:n
                    
                    % The stream function and velocity at this point x
                    psi=-x(n+5)*x(4+i)...
                        +sum(x(n+6:2*n+5)./(js*k).*sinh(k*js.*(x(4+i)+h))...
                        ./cosh(js*k*h).*cos(js*i*pi/n));
                    u=-x(n+5) ...
                        +sum(x(n+6:2*n+5).*cosh(k*js.*(x(4+i)+h))...
                        ./cosh(js*k*h).*cos(js*i*pi/n));
                    w= sum(x(n+6:2*n+5).*sinh(k*js.*(x(4+i)+h))...
                        ./cosh(js*k*h).*sin(js*i*pi/n));
                    
                    % The kinematic condition:  psi(eta)=q
                    F(4+i)=psi-x(2);
                    
                    % The dynamic condition:  g eta + 1/2(u^2+w^2)=R
                    F(4+n+1+i)=g*x(4+i) + 1/2*(u^2 + w^2) - x(3);
                end
                
            otherwise
                error('intype must be either "T" or "L".');
        end
    end
end