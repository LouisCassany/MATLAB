function [sys,x0,str,ts,simStateCompliance] = kite_S_function(t,x,u,flag,x0)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%

switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(x0);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(x0_init)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
sizes = simsizes;

sizes.NumContStates  = 5;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 8;
sizes.NumInputs      = 3;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = x0_init;

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(~,x,u)

[~,params]= init_model();

% Model parameters
r = params(1);
m = params(2);
g_k = params(3);
M_k = params(4);
A_k = params(5);
rho_air = params(6);
g = params(7);
K_i = params(8);
alpha_i0 = params(9);
C_D = params(10);
C_L = params(11);
C_T = params(12);

% Transformation matrices
Mx = @(alpha) [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

% State
theta = x(1);
phi = x(2);
dtheta = x(3);
dphi = x(4);
psi = x(5);

% Input
delta = u(1);
epsilon = u(2);
V_WR = u(3); % disturbance

WR_M_k0 = Mz(phi)*My(theta-pi/2); % from Rk0 to RWR
V_WR = [V_WR;0;0];
V_k_k0 = [-r*dtheta;r*dphi*cos(theta);0]; % Kite's speed in Rk0
V_a_WR = V_WR-WR_M_k0*V_k_k0;
modV_a = norm(V_a_WR);

k0_M_b = Mz(pi+psi);
xb_k0 = k0_M_b*[1;0;0];
yb_k0 = k0_M_b*[0;1;0];
zb_k0 = k0_M_b*[0;0;1];
xa_WR = -V_a_WR/modV_a;
xa_k0 = WR_M_k0.'*xa_WR;

%%

alpha_0 = pi/2 - acos(dot(zb_k0,xa_k0)/(norm(zb_k0)*norm(xa_k0)));
alpha_i = K_i*epsilon+alpha_i0;
alpha = alpha_0 + alpha_i;

% x_ap_k0 = [xa_k0(1);xa_k0(2);0]/norm([xa_k0(1);xa_k0(2);0]);

% C = dot(x_ap_k0,xb_k0)/(norm(dot(x_ap_k0,xb_k0)));
% S = (xb_k0(1)*x_ap_k0(2)-xb_k0(2)*x_ap_k0(1));
% alpha_y = sign(S)*acos(C);

lift = 0.5*rho_air*C_L*alpha*A_k*modV_a*modV_a;
drag = 0.5*rho_air*C_D*alpha*A_k*modV_a*modV_a;
% trans = 0.5*rho_air*C_T*alpha_y*A_k*modV_a*modV_a*1000;


x_lift_k0=cross(-xa_k0,yb_k0)/norm(cross(-xa_k0,yb_k0));
x_drag_k0= -xa_k0;

F_a_k0 = lift*(x_lift_k0)+drag*(x_drag_k0);%+trans*(-yb_k0);

P_k0 = [m*g*cos(theta); 0; m*g*sin(theta)];

F_ext_k0 = P_k0 + F_a_k0;

V_a_Rb=k0_M_b'*WR_M_k0'*V_a_WR;
va = -V_a_Rb(1);

ddtheta = (-F_ext_k0(1)/(r*m))-sin(theta)*cos(theta)*dphi^2;
ddphi = (F_ext_k0(2)/(r*m*cos(theta)))+ 2*tan(theta)*dtheta*dphi;
dpsi = g_k*va*delta+(M_k*((cos(theta)*sin(psi))/va))-dphi*sin(theta);

sys = [dtheta; dphi; ddtheta; ddphi; dpsi];

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

[~,params]= init_model();

% Model parameters
r = params(1);
m = params(2);
g_k = params(3);
M_k = params(4);
A_k = params(5);
rho_air = params(6);
g = params(7);
K_i = params(8);
alpha_i0 = params(9);
C_D = params(10);
C_L = params(11);


% Transformation matrices
Mx = @(alpha) [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

% State
theta = x(1);
phi = x(2);
dtheta = x(3);
dphi = x(4);
psi = x(5);

% Input
delta = u(1);
epsilon = u(2);
V_WR = u(3); % disturbance

WR_M_k0 = Mz(phi)*My(theta-pi/2); % from Rk0 to RWR
V_WR = [V_WR;0;0];
V_k_k0 = [-r*dtheta;r*dphi*cos(theta);0]; % Kite's speed in Rk0
V_a_WR = V_WR-WR_M_k0*V_k_k0;
modV_a = norm(V_a_WR);

k0_M_b = Mz(pi+psi);
xb_k0 = k0_M_b*[1;0;0];
yb_k0 = k0_M_b*[0;1;0];
zb_k0 = k0_M_b*[0;0;1];
xa_WR = -V_a_WR/modV_a;
xa_k0 = WR_M_k0.'*xa_WR;

alpha_0 = pi/2 - acos(dot(zb_k0,xa_k0)/(norm(zb_k0)*norm(xa_k0)));
alpha_i = K_i*epsilon+alpha_i0;
alpha = alpha_0 + alpha_i;

lift = 0.5*rho_air*C_L*alpha*A_k*modV_a*modV_a;
drag = 0.5*rho_air*C_D*alpha*A_k*modV_a*modV_a;

sys = [x; alpha; lift; drag];

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
