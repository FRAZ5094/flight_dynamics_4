%                     ********
%                     * trim *
%                     ********
%
% This programme calculates the longitudinal trim state of business jet
% configuration.  Data is based on the Hawker Siddley HS125 aircraft.  The
% later versions of this were manufactured by Hawker-Beechcraft as the Hawker
% 750/850XP/950 series.
%
% The three longitudinal equations of motion are solved for thrust, pitch
% attitude and elevator angle to fly at a specific speed and altitude
% (given by air density)
%
% Initial guess at solution
%
th=12000;    theta=0.01;   deltae=-0.01;
%
% tol = error tolerence for iteration, step = step size for numerical
% differentiation in Jacobian calculation, pertsmall = smallest allowable 
% perturbaton size (to avoid rounding errors) maxiter = maximum number of
% iterations
%
tol=[0.01 0.01 0.01]; step=0.2; maxiter=10;
pertsmall=[100 0.001 0.001];
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
rho=input('Enter the air density (kg/m^3) > ');
vf=vfk*0.5148;          % Conversion from knots to m/s
%
% Start of iterative loop
%
for iter=1:maxiter
%
% Latest estimate of the unknown variables
%
thi=th;    thetai=theta;   deltaei=deltae; 
%
% The following code perturbs each of the 3 unknowns twice (positive
% pertubation and a negative perturbation) from the current estimate, then 
% calculates the three functions at each perturbation.  The functions are
% also calculated at the unperturbed values (making a total of 7 sets of 
% calculations. Central differencing is then used to differentiate the 
% functions w.r.t. the unknowns to construct the jacobian.
%
var(1)=thi;
var(2)=thetai;
var(3)=deltaei;
%
% Calculate perturbation sizes
%
for lv=1:3
    pert(lv)=abs(var(lv))*step;
    if pert(lv)<pertsmall(lv) pert(lv)=pertsmall(lv); end
    dvar(lv)=2*pert(lv);
end
%
% Set up the values of the 3 unknowns for the 7 calculations
%
for ncalc=1:7
    for nvar=1:3
        varpert(ncalc,nvar)=var(nvar);
    if ncalc>1 
        np=fix(0.5*ncalc);
        if rem(ncalc,2)==1 
            varpert(ncalc,np)=var(np)-pert(np);
        else
            varpert(ncalc,np)=var(np)+pert(np);
        end
    end
    end
end
%
% Now calculate the functions
%
for j=1:7
%
% "Unknown" variables
%
    th=varpert(j,1);
    theta=varpert(j,2);
    deltae=varpert(j,3);
%
% Aircraft States
%
u=vf*cos(theta);
w=vf*sin(theta);
q=0;                 % No angular velocity in trim
%
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
% Lift, Darag, Pitching Moment
%
l=const*cl;
d=const*cd;
ma=const*cbar*cm;
%
% External Forces and Moments
%
x=l*sin(alpha)-d*cos(alpha)+th;
z=-l*cos(alpha)-d*sin(alpha);
mp=ma+th*ht;
%
% Functions (equations of motion)
%
f(1,j)=x-mass*9.81*sin(theta);
f(2,j)=z+mass*9.81*cos(theta);
f(3,j)=mp;
end
%
% Calculate the Jacobian by central differences
%
for mfunc=1:3
    for mvar=1:3
        jac(mfunc,mvar)=(f(mfunc,mvar*2)-f(mfunc,(2*mvar+1)))/dvar(mvar);
    end
end
%
% Invert the Jacobian
%
jaci=inv(jac);
%
jif=jaci*f(:,1);
%
deltae=deltae+pert(3);
%
% Calculate new estimates of the unknowns
%
th=thi-jif(1);
theta=thetai-jif(2);
deltae=deltaei-jif(3);
%
% Calulate the error and check for convergence
%
err(1)=abs((th-thi)/th);
err(2)=abs((theta-thetai)/theta);
err(3)=abs((deltae-deltaei)/deltae);
%
if err<tol break
end
%
% end of iterative loop
%
end
%
% Output to screen
%
fprintf('\n')
fprintf('Thrust = %.3f kN\n',th/1000)
fprintf('Pitch Attitude = %.4f rad  (%.2f deg)\n',theta,theta/0.01745)
fprintf('Elevator Angle = %.4f rad  (%.2f deg)\n',deltae,deltae/0.01745)
