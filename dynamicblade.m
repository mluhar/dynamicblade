% This MATLAB script numerically simulates the dynamics of a flexible blade
% under oscillatory flow.  The blade is modeled as an inextensible, linear
% elastic beam undergoing finite (i.e., nonlinear deformation).  The fluid
% forces are modeled based on the well-known Morrison force formulation.

% This is a finite difference solver that is 2nd order accurate in space.
% The blade-normal force balance is solved explicitly to yield the bending 
% angle, theta, while the blade-tangential force balance is solved 
% implicitly to yield tension.  

% For further details, see:
% 1. M Luhar and HM Nepf (2016) Wave-induced dynamics of flexible blades,
% Journal of Fluids and Structures, 61: 20-41
% 2. M Luhar (2012) Analytical and experimental studies of plant-flow
% interaction at multiple scales, PhD thesis, MIT

% The publications above used PIV-measured velocity fields to force the
% blade, here we use a sinusoidal velocity field. 

% Developed by M Luhar* and HM Nepf+
% * Department of Aerospace and Mechanical Engineering, USC
% + Department of Civil and Environmental Engineering, MIT

% Written by M Luhar(luhar@usc.edu, mluhar@alum.mit.edu, mluhar@cantab.net)
% Please contact me in case of any questions or concerns. 

% SEE ATTACHED LICENSING AGREEMENT AND README FILE

%% Dimensionless Inputs (See Luhar and Nepf 2016)
rhop = 1;       %Density ratio between blade and water
B = 1;          %Buoyancy Parameter: Ratio of Buoyancy to Rigidity
Ca = 1e3;       %Cauchy Number: Ratio of Hydrodynamic Drag to Rigidity
KC = 20;        %Keulegan Carpenter Number: Ratio of Inertia to Drag
L = 4;          %Ratio of blade length to orbital excursion
S = 0.01;       %Ratio of blade thickness to width
CM = 1;         %Added Mass Coefficient
CD = max(1.95,10*KC^(-1/3)); %Drag Coefficient
CF = 0.2;       %Friction Coefficient

%% Set up model
% Model parameters
% NOTE: You will probably have to tinker with the number of grid points and
% time step to ensure stability.
ns = 256;               %Number of grid points over s in [0,1]
s  = linspace(0,1,ns)'; %Grid: s represents distance along blade from base
ds = s(2)-s(1);         %spacing
dt = 0.10;              %Time
nt = ceil(10*2*pi/dt);  %Number of time steps for simulation
i = sqrt(-1);

%2nd order accurate difference matrices. 
D1 = fdmatrix(s,1,2);
D2 = fdmatrix(s,2,2);
D3 = fdmatrix(s,3,2);

%Make matrices for finite differences
%Make matrix for implicit Tension calculation
D1T = D1;
%Fixed BC
D1T(end,:) = 0; D1T(end,end) = 1;  
%Make matrix for conversion between theta and X
D1X = D1;
%Fixed BC
D1X(1,:) = 0; D1X(1,1) = 1;     

%% Initial conditions
%velocity and acceleration
u = zeros(ns,1);
ut = zeros(ns,1);

%theta
theta = zeros(ns,1);
thOld = theta;          %Previous time step
thetat = zeros(ns,1);   %time derivative

%Calculate X from theta
RHSX = i*exp(-i*theta); RHSX(1)=0;
X = D1X\RHSX;
XOld = X;               %Previous time step
XOld2= X;               %Two time steps ago
Xt = zeros(ns,1);       %time derivative

%Tension
T = linspace(B,0,ns)';

%% March in time
for tc = 1:nt
    %Water Velocity and Acceleration, time is normalized by omega
    %Assume sinusoidal for now
    u = (1-exp(-tc*dt))*sin(tc*dt);
    ut= (1-exp(-tc*dt))*cos(tc*dt);
    
    %Relative velocity normal to blade
    UN = abs(real(exp(i*theta).*(u-L*Xt)));
    
    %Make lower triangular matrices for integration
    IS = ds*tril(sin(theta*ones(1,ns)-(theta*ones(1,ns))'),-1);
    IC = ds*tril(cos(theta*ones(1,ns)-(theta*ones(1,ns))'),-1);
    
    %Calculate all the terms treated explicitly
    %Tension, Buoyancy
    A1 = T.*(D1*theta) - B*sin(theta);
    %Drag Forcing
    A2 = (1/2)*Ca*CD.*UN.*real(exp(i*theta).*u);
    %Added mass and virtual buoyancy forcing
    A3 = Ca*real(exp(i*theta).*((2*pi/KC)*(pi*CM/4+S).*ut));
    %Terms due to time discretization
    A4 = -Ca*L*(2*pi/KC)*(pi*CM/4+rhop*S).*(IS*(thetat.^2));
    A5 = (1/2)*Ca*L*CD.*UN.*(IC*(thOld/2/dt));
    A6 = Ca*L*(IC*((2*pi/KC)*(pi*CM/4+rhop*S).*(2*theta-thOld)/dt/dt));
    
    %Total Forcing
    RHS = (A1+A2+A3+A4+A5+A6);
    %Account for BCs
    RHS(1) = 0;
    RHS(end) = 0;
    RHS(end-1) = 0;
    
    %Calculate the matrix for the LHS
    LHS = (D3 + Ca*L*diag((1/2)*CD.*UN/2/dt + (2*pi/KC)*(pi*CM/4+rhop*S)/dt/dt)*IC);
    %Fixed boundary condition
    LHS(1,:) = 0;
    LHS(1,1) = 1; 
    %No bending moment
    LHS(end,:) = D1(end,:); 
    %No shear force
    LHS(end-1,:) = D2(end,:); 
    
    %Evaluate!
    thNew = LHS\RHS;
    RHSX = i*exp(-i*thNew);
    RHSX(1)=0;
    XNew = D1X\RHSX;
    
    %Update variables
    %Rates of change
    Xt = (1.5*XNew-2*X+0.5*XOld)/dt;
    Xtt = (2.0*XNew - 5*X + 4*XOld - XOld2)/dt/dt;
    thetat = (1.5*thNew-2*theta+0.5*thOld)/dt;
    %Theta
    thOld = theta;
    theta = thNew;
    %X
    XOld2= XOld;
    XOld = X;
    X = XNew;
    
    %Relative velocity along blade
    UT = abs(imag(exp(i*theta).*(u-L*Xt)));

    %Now calculate tension implicitly;
    B1 = -(D2*thNew).*(D1*thNew);
    B2 = -B*cos(thNew);
    B3 = -Ca*(2*pi/KC)*S*imag(exp(i*thNew).*(ut-rhop*L*Xtt));
    B4 = -(1/2)*CF*Ca*UT.*imag(exp(i*thNew).*(u-L*Xt));
    RHST = B1+B2+B3+B4;
    RHST(end) = 0;
    T = D1T\RHST;

    %Plot new solution every 5 time steps
    if(mod(tc,5)==0)
        figure(1)
        clf
        plot(real(X),imag(X),'ro-','linewidth',1,'markersize',3);
        hold on
        quiver(-0.5,1.25,real(u(end)),imag(u(end)),0.5,'k','linewidth',2)
        pbaspect([2 1.5 1]); xlim([-1 1]); ylim([0 1.5]);
        title(strcat('step=',num2str(tc),', time=',num2str(tc*dt)),'fontsize',16); drawnow
    end

    %Parameters to be saved
    saved.t(:,tc) = tc*dt;      % time
    saved.theta(:,tc) = theta;  % bending angle
    saved.X(:,tc) = X;          % position
    saved.U(:,tc) = u;          % velocity
    saved.T(:,tc) = T;          % tension
    D2th = D2*theta;            
    saved.F(tc) = D2th(1);      % force
end