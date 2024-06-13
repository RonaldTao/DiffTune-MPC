%
% Copyright (c) The acados authors.
%
% This file is part of acados.
%
% The 2-Clause BSD License
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.;

%

function model = unicycle_model()

import casadi.*

%% dims
nx = 3; % four states: px, py, vx, vy
nu = 2; % two controls: ax,ay

Ts = 0.05;
%% symbolic variables
sym_x1 = SX.sym('x', 1, 1); % x position
sym_y = SX.sym('y', 1, 1); % y position
sym_theta = SX.sym('theta', 1, 1); % angle
sym_x = vertcat(sym_x1,sym_y,sym_theta); % states

sym_u = SX.sym('u', nu, 1); % controls

sym_x1dot = SX.sym('x1dot', 1, 1); % x position
sym_ydot = SX.sym('ydot', 1, 1); % x position
sym_thetadot = SX.sym('thetadot', 1, 1); % x position

sym_xdot = vertcat(sym_x1dot,sym_ydot,sym_thetadot); % states



%% dynamics
x1_dynamics = cos(sym_theta)*sym_u(1); %continous   
x1_kp1 = sym_x1+x1_dynamics*Ts;     % discrete

y_dynamics = sin(sym_theta)*sym_u(1);
y_kp1 = sym_y+y_dynamics*Ts;

theta_dynamics = sym_u(2);
theta_kp1 = sym_theta+theta_dynamics*Ts;

expr_f_expl=vertcat(x1_dynamics,y_dynamics,theta_dynamics);

expr_f_impl = expr_f_expl-sym_xdot;

expr_f_discrete = vertcat(x1_kp1,y_kp1,theta_kp1);

%% constraints
expr_h = [sym_u; sym_x];
expr_h_e = [sym_x];

%% nonlnear least squares
expr_y = [sym_u; sym_x];
expr_y_e = [sym_x];

%% external cost (Sheng: not clear what this is doing.)
% yr_u = zeros(nu, 1);
% yr_x = zeros(nx, 1);
% dWu = 2*ones(nu, 1);
% dWx = ones(nx, 1);
% 
% ymyr = [sym_u; sym_x] - [yr_u; yr_x];
% ymyr_e = sym_x - yr_x;
% 
% expr_ext_cost = 0.5 * ymyr' * ([dWu; dWx] .* ymyr);
% expr_ext_cost_e = 0.5 * ymyr_e' * (dWx .* ymyr_e);

%% populate structure
d = 0.5; % meter, distance between two wheels
r = 0.1; % meter, radius of each wheel
w_min = -2*pi;
w_max = 2*pi;
model.d = d;
model.r = r;
model.w_min = w_min;
model.w_max = w_max;
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;
%model.expr_phi = expr_phi;
model.expr_h = expr_h;
model.expr_h_e = expr_h_e;
model.Ts = Ts; % add sample time
model.expr_f_discrete = expr_f_discrete;
%model.A = A; % save the discrete-time system matrices
%model.B = B; % save the discrete-time system matrices
% model.expr_y = expr_y;
% model.expr_y_e = expr_y_e;
% model.expr_ext_cost = expr_ext_cost;
% model.expr_ext_cost_e = expr_ext_cost_e;
