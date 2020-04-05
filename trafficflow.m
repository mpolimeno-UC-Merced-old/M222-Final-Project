%code to model traffic flow

%traffic at  stop light
%modifying code from garcia II edition

clear all;
N = 80;
L = 400; %meters of track

h = L/N;

vmax = 25;%maximum car speed

k = 0.2;%h/vmax; %to simulate realistic scenario

%set moment in time when last car starts moving

fprintf('Last car starts moving after %g steps\n',(L/4)/(vmax*k));
nstep = input('Enter number of steps: '); %last car starts moving when traffic is no longer bumper to bumper

mu = k/(2*h); %for FTCS part
lambda = k^2/(2*h^2); %for Lax-Wendroff

%set IC's and periodic boundary
%periodic boundaries are chosen to resemble the start of a race on a
%circular track
umax = 1; %max density
Qmax = 0.25*umax*vmax; %flux=density x speed

u = zeros(1,N);
for ii=round(N/4):round(N/2-1)
    u(ii) = umax;
end

%periodic BC's
ip(1:N) = (1:N)+1; %position i+1 at time n
im(1:N) = (1:N)-1; %postion i-1 at time n
ip(N) = 1;
im(1) = N;
jj = 1;
while jj<=nstep
    Q = u.*(vmax*(1-u/umax));
    %lax-wendroff
    cpstar = vmax*(1- (u(ip)+u(1:N))/umax);
    cmstar = vmax*(1- (u(1:N)+u(im))/umax);
    u(1:N) = u(1:N) - mu*(Q(ip) - Q(im)) + lambda*(cpstar.*(Q(ip)+Q(1:N))...
        -cmstar.*(Q(1:N)-Q(im)));
    jj = jj + 1;
end

Nvec = 1:N;
plot(Nvec,u);