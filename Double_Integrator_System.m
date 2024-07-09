% This script defines the system model for a double integrator system.
% More info about the system can be found in Section VI of the paper
% 'DiffTune-MPC: Closed-Loop Learning for Model Predictive Control' (https://arxiv.org/pdf/2312.11384)

% Ran Tao, Sheng Cheng
% University of Illinois Urbana-Champaign

function model = Double_Integrator_System()

import casadi.*

%% dims
nx = 2;
nu = 1;

%% symbolic variables
sym_x = SX.sym('x', nx, 1); % states
sym_u = SX.sym('u', nu, 1); % controls
sym_xdot = SX.sym('xdot',size(sym_x)); %state derivatives

%% dynamics
% continuous time
Ac = [0 1;0 -0.05];
% 0.05 is the friction efficient, assumed to be proportional to velocity
% but in reverse direction

Bc = [0;1];

% discrete time
Ts = 0.01; % sampling time
M = expm([Ts*Ac, Ts*Bc; zeros(nu, 2*nx/2+nu)]);
A = M(1:nx,1:nx);
B = M(1:nx,nx+1:end);

expr_f_expl = Ac*sym_x + Bc*sym_u;
expr_f_impl = expr_f_expl - sym_xdot;
expr_phi = A*sym_x + B*sym_u;

u_max = 1; %maximum control input
x_max = 100; %maximum position
v_max = 100; %maximum velocity
%% constraints
expr_h = [sym_u; sym_x];
expr_h_e = [sym_x];

%% nonlnear least squares
expr_y = [sym_u; sym_x];
expr_y_e = [sym_x];

%% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_phi = expr_phi;
model.expr_h = expr_h;
model.expr_h_e = expr_h_e;
model.Ts = Ts; % add sample time
model.A = A; % save the discrete-time system matrices
model.B = B; % save the discrete-time system matrices
model.u_max = u_max;
model.x_max = x_max;
model.v_max = v_max;
