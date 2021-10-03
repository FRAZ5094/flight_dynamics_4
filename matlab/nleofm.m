function [ydot]=nleofm(t,y)

global mass iyy s cbar ht cd0 cda cdaa cl0 cla clde cm0 cma cmde cmq rho... 
    deltaee deltaestep the thestep;
%
% Assign latest values of y to the states and controls
%
u=y(1);
w=y(2);
q=y(3);
theta=y(4);
%
deltae=deltaee+deltaestep;  % Add step in elevator to trim value
thrust=the+thestep;   % Add thrust step to trim value

vf=sqrt(u*u+w*w);

alpha=atan(w/u);
qh=q*cbar/vf;
%
% Aerodynamic coefficients
%
cl=cl0+cla*alpha+clde*deltae;
cd=cd0+cda*alpha+cdaa*alpha*alpha;
cm=cm0+cma*alpha+cmde*deltae+cmq*qh;
%
const=0.5*rho*vf*vf*s;
%
% Lift, Drag, Pitching Moment
%
l=const*cl;
d=const*cd;
ma=const*cbar*cm;
%
% External Forces and Moments
%
x=l*sin(alpha)-d*cos(alpha)+thrust;
z=-l*cos(alpha)-d*sin(alpha);
mp=ma+thrust*ht; 
%
% Equations of motion in standard (ydot = f(y)) form
%
ydot=zeros(6,1);    % Ensure ydot is recognised as a column matrix
%
ydot(1)=-q*w+x/mass-9.81*sin(theta);
ydot(2)=q*u+z/mass+9.81*cos(theta);
ydot(3)=mp/iyy;
ydot(4)=q;
ydot(5)=u*cos(theta)-w*sin(theta);
ydot(6)=u*sin(theta)+w*cos(theta);
