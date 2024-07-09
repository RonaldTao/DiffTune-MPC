% This script defines the system model for differential wheeled robot

% Ran Tao, Sheng Cheng
% University of Illinois Urbana-Champaign

function model = Differential_Wheeled_Robot()

import casadi.*

%% dims
nx = 3; % four states: px, py, vx, vy
nu = 2; % two controls: ax,ay
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
ng = 1; %when we have linear inequality constraints on state and input, set ng = 1
nh = 0;
%% symbolic variables
sym_x1 = SX.sym('x', 1, 1); % horizontal position
sym_y = SX.sym('y', 1, 1); % vertical position
sym_theta = SX.sym('theta', 1, 1); % angle
sym_x = vertcat(sym_x1,sym_y,sym_theta); % total state
sym_u = SX.sym('u', nu, 1); % control input
sym_x1dot = SX.sym('x1dot', 1, 1); % horizontal position time derivative
sym_ydot = SX.sym('ydot', 1, 1); % vertical position time derivative
sym_thetadot = SX.sym('thetadot', 1, 1); % angle time derivative
sym_xdot = vertcat(sym_x1dot,sym_ydot,sym_thetadot); % total state derivative

%% dynamics
x1_dynamics = cos(sym_theta)*sym_u(1); %continous time dynamics
Ts = 0.05; %% System sampling time for discrete system
x1_kp1 = sym_x1+x1_dynamics*Ts;     % discrete time dynamics
y_dynamics = sin(sym_theta)*sym_u(1);
y_kp1 = sym_y+y_dynamics*Ts;
theta_dynamics = sym_u(2);
theta_kp1 = sym_theta+theta_dynamics*Ts;
expr_f_expl = vertcat(x1_dynamics,y_dynamics,theta_dynamics); %explicit expression for continous time
expr_f_impl = expr_f_expl - sym_xdot; %implicit expression for continous time
expr_f_discrete = vertcat(x1_kp1,y_kp1,theta_kp1); %discrete time

%% constraints
expr_h = [sym_u; sym_x];
expr_h_e = [sym_x];

%% nonlnear least squares
expr_y = [sym_u; sym_x];
expr_y_e = [sym_x];

%% populate structure
d = 0.5; % meter, distance between two wheels
r = 0.1; % meter, radius of each wheel
w_min = -2*pi; %minimum angular velocity of each wheel
w_max = 2*pi; %maximum angular velocity of each wheel

model.d = d;
model.r = r;
model.w_min = w_min;
model.w_max = w_max;
model.nx = nx;
model.nu = nu;
model.ny = ny;
model.ny_e = ny_e;
model.ng = ng;
model.nh = nh;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
model.expr_h = expr_h;
model.expr_h_e = expr_h_e;
model.Ts = Ts;
model.expr_f_discrete = expr_f_discrete;
