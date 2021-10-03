%              ***********
%              * bjetsim *
%              ***********
%
% This program performs a time response calculation for an elevator step
% input to a simple nonlinear longitudinal model of a business jet 
% configuration.  Data is based on the Hawker Siddley HS125 aircraft.  The
% later versions of this were manufactured by Hawker-Beechcraft as the Hawker
% 750/850XP/950 series.
%
% Have to make data available in functions
%
global mass iyy s cbar ht cd0 cda cdaa cl0 cla clde cm0 cma cmde cmq rho... 
    deltaee deltaestep the thestep;
%
% Basic Aircraft Data
%
mass=7485.0; iyy=84309;
s=32.8;      cbar=2.29;   ht=-0.387;   
cd0=0.177;   cda=0.232;   cdaa=1.393;
cl0=0.895;   cla=5.01;    clde=0.722;
cm0=-0.046;  cma=-1.087;  cmde=-1.88;  cmq=-7.055;
%
% User Input
%
vfk=input('Enter the flight speed of the aircraft (knots)  > ');
vf=vfk*0.5148;          % Conversion from knots to m/s
rho=input('Enter the air density (kg/m^3) > ');
thetae=input('Enter the trim pitch attitude (radians)  > ');
the=input('Enter the thrust in the trim state (N) > ');
thestep=input('Enter the amplitude of the thrust step input (N)  >');
deltaee=input('Enter the elevator angle for trim (radians) > ');
deltaestep=input('Enter the amplitude of the elevator step input (radians)  > ');
tmax=input('Enter the total time for simulation > ');
fname=input('Enter file name for output >  ','s');
%
% Calculate initial conditions for the integration (i.e. the trim state)
%
ue=vf*cos(thetae);
we=vf*sin(thetae);
qe=0; xe=0; ze=0;
%
y0=[ue, we, qe, thetae, xe, ze];
%
% Solve the equations of motion using a Runge-Kutta routine
%
[t,y]=ode45('nleofm',[0 tmax],y0);
%
% Write Data to File
%
save (fname, 't', 'y')
%
% Plot responses
%
subplot(3,2,1)
plot(t,y(:,1))
xlabel('Time (s)')
ylabel('U(m/s)')
subplot(3,2,2)
plot(t,y(:,2))
xlabel('Time (s)')
ylabel('W(m/s)')
subplot(3,2,3)
plot(t,y(:,3)/0.01745)
xlabel('Time (s)')
ylabel('Q(deg/s)')
subplot(3,2,4)
plot(t,y(:,4)/0.01745)
xlabel('Time (s)')
ylabel('Theta(deg)')
subplot(3,2,5)
plot(y(:,5),y(:,6))
xlabel('x (m)')
ylabel('Altitude (m)')


