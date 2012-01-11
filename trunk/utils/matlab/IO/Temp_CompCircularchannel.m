close all
%
% This script is useful for visualizing data output (in binary file format)
% from the fortran SWENSE3D model. The script can be easily tailored
% to output data from other data files as well.
%
% By Allan P. Engsig-Karup, apek@imm.dtu.dk.
clear all
% Show free surface evolution
close all

curdir = cd;

% Path for data files here
dirpath = '..';
cd(dirpath)

initialstep = 1;
Nsteps = 641; %1001; 
jump   = 1;
dt     = 0.0234; %1.09414167790990047E-2;%1.08328305997448821E-2; %1.171537321091214E-002;% 1;%1.08328305997448821E-2;
g      = 9.81;
plotmethod = 3;
IOmethod = 2; %0:binary ; 1:classical unformatted ; 2:unformatted ftn95
fac = 10;
i=1;

% CHOOSE BETWEEN byte storage format below (uncomment the one needed)
% byteorder = 'ieee-be'; % IEEE Big-Endian format 
byteorder = 'ieee-le'; % IEEE Little-Endian format

for tstep = initialstep : jump : Nsteps 
    tid = tstep*dt+dt;

    if mod(tstep,10)==0
        disp(sprintf('tstep=%d...',tstep))
    end

    switch IOmethod
        case 1 % unformatted file
            fid=fopen(sprintf('EP_%.5d.bin',tstep),'r',byteorder);
            [nrec,count] = fread(fid,1,'int32');
            Nx=fread(fid,1,'int32');
            Ny=fread(fid,1,'int32');
            [nrec,count] = fread(fid,1,'int32');

            [nrec,count] = fread(fid,1,'int32');
            X=fread(fid,[Nx Ny],'float64');
            Y=fread(fid,[Nx Ny],'float64');
            [nrec,count] = fread(fid,1,'int32');
            
            [nrec,count] = fread(fid,1,'int32');
            E=fread(fid,[Nx Ny],'float64');
            %                         E=reshape(E(:),nx,ny);
            P=fread(fid,[Nx,Ny],'float64');
            fclose(fid);
         case 2 % unformatted file from ftn95
            fid=fopen(sprintf('EP_%.5d.bin',tstep),'r',byteorder);
            [nrec,count] = fread(fid,1,'int8');
            Nx=fread(fid,1,'int32');
            Ny=fread(fid,1,'int32');
            [nrec,count] = fread(fid,1,'int8');
            if(2*Nx*Ny > 30)
                [nrec,count] = fread(fid,1,'int8');
                [nrec,count] = fread(fid,1,'int32');
                X=fread(fid,[Nx Ny],'float64');
                Y=fread(fid,[Nx Ny],'float64');
                [nrec,count] = fread(fid,1,'int32');
                [nrec,count] = fread(fid,1,'int8');

                [nrec,count] = fread(fid,1,'int8');
                [nrec,count] = fread(fid,1,'int32');
                E=fread(fid,[Nx Ny],'float64');
                %                         E=reshape(E(:),nx,ny);
                P=fread(fid,[Nx,Ny],'float64');
                fclose(fid);  
            else
                [nrec,count] = fread(fid,1,'int8');
                X=fread(fid,[Nx Ny],'float64');
                Y=fread(fid,[Nx Ny],'float64');
                [nrec,count] = fread(fid,1,'int8');

                [nrec,count] = fread(fid,1,'int8');
                E=fread(fid,[Nx Ny],'float64');
                %                         E=reshape(E(:),nx,ny);
                P=fread(fid,[Nx,Ny],'float64');
                fclose(fid);  
            end  
        otherwise  % binary file
            fid=fopen(sprintf('EP_%.5d.bin',tstep),'r',byteorder);
            Nx=fread(fid,1,'int32');
            Ny=fread(fid,1,'int32');
            X=fread(fid,[Nx Ny],'float64');
            Y=fread(fid,[Nx Ny],'float64');
            E=fread(fid,[Nx Ny],'float64');
            %                         E=reshape(E(:),nx,ny);
            P=fread(fid,[Nx,Ny],'float64');
            fclose(fid);
    end

% time for plotting...
if Nx>1 & Ny>1
    plotmethod = 5; %2;
%   Use the following for cylinder tests
%     Esave(1,i)=E(2,2); % (\beta = 0)
%     if mod(Nx-2,2) == 0 %even
%         Esave(2,i)=(E(1+(Nx-2)/2,2)+E(1+(Nx-2)/2+1,2))/2; % (\beta = \pi / 2)
%     else %odd
%         Esave(2,i)=E(1+ceil((Nx-2)/2),2); % (\beta = \pi / 2)
%     end
%     Esave(3,i)=E(Nx-1,2); % (\beta = \pi)
%     % nondimensionalize: 2*\eta / H
%     Esave(:,i) = 2*Esave(:,i)/0.01;
%     % time scale...
     time(i) = (i-1)*dt;
     i=i+1;
else
    plotmethod = 1;
    Esave(:,i)=E;
    i=i+1;
end
switch plotmethod
    case 1
        if Nx>1
            plot(X,E,'b')%,X,P,'r')
            axis([X(1) X(end) -0.05 0.1])
        else
            plot(Y(2:Ny-1),E(2:Ny-1),'b',Y(2:Ny-1),P(2:Ny-1),'r')
            %axis([Y(1) Y(end) -0.05 0.1])
            %pause
        end
        grid on
    case 2
        %             mesh(X,Y,Eerr*fac)
        %surf(X,Y,E*fac)
        surf(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),E(2:end-1,2:end-1)*fac)
        colormap(water)
        %axis([X(1) X(end) Y(1) Y(end) -0.1 0.1])
        axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y)) -0.1 0.1])
        light
        %             axis equal
        %axis off
        %             view(0,90)
        % view(-10,0)
    case 3
        subplot(1,3,1)
        mesh(X,Y,E*fac)
        %             axis([x(1) x(end) y(1) y(end) -1 1])
        subplot(1,3,2)
        mesh(X,Y,P*fac)
        %             axis([x(1) x(end) y(1) y(end) -1 1])
        subplot(1,3,3)
        mesh(X,Y,abs(E-P)*fac)
    case 4
        %             mesh(X,Y,Eerr*fac)
        surf(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),(E(2:end-1,2:end-1)+0*P(2:end-1,2:end-1))*fac)
        colormap(water)
        axis([X(1) X(end) Y(1) Y(end) -0.1 0.1])
        light
        %             axis equal
        %axis off
        %             view(0,90)
        % view(-10,0)
     case 5    
        %contourf(X,Y,E)
        v=[-0.1:0.0075:0.1]; % warning, no contour on 0 !!
        contourf(X(2:Nx-1,2:Ny-1),Y(2:Nx-1,2:Ny-1),E(2:Nx-1,2:Ny-1),v)
        caxis([-0.05 0.05])
end
drawnow
% pause(0.3)
end
if Nx>1 & Ny>1
%   Use this for cylinder problem
%     figure(2)
%     subplot(3,1,1)
%     plot(time,Esave(1,:))
%     axis([0 time(end) -2 2])
%     subplot(3,1,2)
%     plot(time,Esave(2,:))
%     axis([0 time(end) -2 2])
%     subplot(3,1,3)
%     plot(time,Esave(3,:))
%     axis([0 time(end) -2 2])
else
    figure(2)
    plot(X(2:Nx-1),Esave(2:Nx-1,size(Esave,2)-50:size(Esave,2)))%,X,P,'r')
    axis([X(2) X(Nx-1) -0.05 0.1])
    figure(3)
    plot(Esave(Nx-1,:))
end
%
% For comparison in the case of semi circular channel
%
cd 'C:\Documents and Settings\guduc\My Documents\Work\3D\SWENSE3D_New\Matlab';
Npts=71;
[x,y,ele, r1, r2, l] = DalrympleEtAll_SemiCircularChannelGeneral(0.5,2,1,Npts);
i=1;
% WAVE EVOLUTION
for t = 0:0.05:10
    %surf(x,y,real(ele*exp(-sqrt(-1)*t)))
    %axis([0 r2 -r2 r2 -max(abs(ele(:))) max(abs(ele(:)))])
    %view(-13,76)
    %drawnow
    ref(1,i)=real(ele(1,1)*exp(-sqrt(-1)*t));% (\beta = 0)
    if mod(Npts,2) == 0 %even
        ref(2,i)=real((ele(1+(Npts)/2,1)+ele(1+(Npts)/2+1,1))*exp(-sqrt(-1)*t));  % (\beta = \pi / 2)
    else %odd
        ref(2,i)=real((ele(1+ceil(Npts/2),1))*exp(-sqrt(-1)*t)); % (\beta = \pi / 2)
    end
    ref(3,i)=real(ele(Npts,1)*exp(-sqrt(-1)*t));% (\beta = \pi)
    time2(i)=t/sqrt(g*2*pi/l);
    i=i+1;
end
figure(3)
subplot(3,1,1)
plot(time2,ref(1,:))
axis([0 time2(end) -2 2])
subplot(3,1,2)
plot(time2,ref(2,:))
axis([0 time2(end) -2 2])
subplot(3,1,3)
plot(time2,ref(3,:))
axis([0 time2(end) -2 2])

figure(4) % warning not done in the same propagation direction...
plot(-y(Npts,:),real(ele(Npts,:)*exp(-sqrt(-1)*(sqrt(g*2*pi/l)*time(end)))))
hold on
plot(X(:,2),2*E(:,2)/0.01,'r--')
hold off

figure(5)% warning not done in the same propagation direction...
plot(-y(1,:),real(ele(1,:)*exp(-sqrt(-1)*(sqrt(g*2*pi/l)*time(end)))))
hold on
plot(X(:,Ny-1),2*E(:,Ny-1)/0.01,'r--')
hold off
