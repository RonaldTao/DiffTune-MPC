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

%% example of closed loop simulation
clear all
clc
close all
% This script implements Difftune-MPC on a differential wheeled robot.
% The original MPC problem is solved using acaods and the analytical gradients are obtained by solving extra MPC problems (LMPC-Grad) using quadprog since the linearized system is time-varying.
% Please install acados in matlab before running this example.
% We define constraints, cost functions and everything in this example following instructions from acados.
% For more info, check problem_formulation_ocp_mex.pdf from acados

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    ! source env.sh
    % 	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

%% setup problem
model = Differential_Wheeled_Robot; %achieve model

% dims
T = 10.0; % total horizon time for the simulation
nx = model.nx; % number of states
nu = model.nu; % number of inputs
ny = model.ny; % number of outputs in lagrange term
ny_e = model.ny_e; % number of outputs in mayer term
ng = model.ng;
nh = model.nh;
% constraint formulation
% bounds on x and u


%% Set up acados
ocp_N = 10; % prediction horizon for discrete time system in MPC
% Initial weight matrix in lagrange term
W = diag([1*ones(nx,1);... % pos, vel
            1*ones(nu,1) ]); % control
[ocp,sim] = setup_acados(model,ocp_N,W);

%% add desired trajectory for tracking
pos_1_des = @(t) 1-cos(t/2);
pos_2_des = @(t) 0.5*t;
pos_3_des = @(t) atan2(0.5,0.5*sin(t/2));

%% Difftune-MPC
% Parameters for DiffTune
learningRate = 0.01; 
minCoef = 0.01; % the minimial diagonal elements of Q and R cannot be smaller than this value
maxCoef = 1000; % the minimial diagonal elements of Q and R cannot be bigger than this value

loss_hist = []; % loss history over iteration
param_hist = [diag(W)]; %parameter history over iteration
RMSE_hist = []; %RMSE history over iteration
theta_hist = []; %parameter theta history over iteration
u_hist = []; %control input history over iteration

% initialize the gradient
theta_gradient = zeros(1,nx+nu); 
theta_gradient_hist = theta_gradient;
total_itr=20; %total iteration for DiffTune
difftune_itr = 0; %initialize the iteration index to 0

%calculate the jacobian of the system dynamics w.r.t. state and control input
import casadi.*
grad_f_x=jacobian(model.expr_f_discrete,model.sym_x);
grad_f_x_fcn = Function('grad_f_x_fcn',{model.sym_x,model.sym_u},{grad_f_x});

grad_f_u=jacobian(model.expr_f_discrete,model.sym_u);
grad_f_u_fcn = Function('grad_f_u_fcn',{model.sym_x,model.sym_u},{grad_f_u});

%% loop
while(1)
    theta_gradient = zeros(1,nx+nu);
    difftune_itr = difftune_itr + 1;
    %% closed loop simulation
    n_sim = T/model.Ts;
    x_sim = zeros(nx, n_sim+1);
    x_sim(:,1) = [0;0;pi/2]; %initial state
    u_sim = zeros(nu, n_sim);
    
    % set the initial guess of x, u, and pi (costate)
    x_traj_init = zeros(nx, ocp_N+1);
    u_traj_init = zeros(nu, ocp_N);
    pi_traj_init = zeros(nx, ocp_N);
    
    desired_state = [pos_1_des(model.Ts * [0:n_sim+ocp_N-1]);
                    pos_2_des(model.Ts * [0:n_sim+ocp_N-1]);
                    pos_3_des(model.Ts * [0:n_sim+ocp_N-1])];

    % initialize the sensitivity
    dx_dtheta = zeros(nx,nx+nu);

    % for each iteration in the following loop, we solve the MPC problem using acados, and achieve the
    % dx/dtheta and du/dtheta in equation (5)
    % After the loop, we have run simulation over entire time hoziron T, and then use
    % the calculated dx/dtheta and du/dtheta to update theta, which is Q
    % and R in W
    
    for ii=1:n_sim
        % desired_state(index,time)
        ocp.set('constr_x0', x_sim(:,ii));
        ocp.set('init_x', x_traj_init);
        ocp.set('init_u', u_traj_init);
        ocp.set('init_pi', pi_traj_init);
        
        % set the y_ref and y_ref_e for tracking a given signal
        time_varying_y_ref = [pos_1_des((ii+[0:ocp_N-1])*model.Ts);
                              pos_2_des((ii+[0:ocp_N-1])*model.Ts);
                              pos_3_des((ii+[0:ocp_N-1])*model.Ts);
                              zeros(nu,ocp_N)];
        for jj = 1:ocp_N
            ocp.set('cost_y_ref', time_varying_y_ref(:,jj), jj-1);
        end
        time_varying_y_ref_N = [pos_1_des((ii+ocp_N)*model.Ts);
                            pos_2_des((ii+ocp_N)*model.Ts);
                            pos_3_des((ii+ocp_N)*model.Ts)];
        ocp.set('cost_y_ref_e', time_varying_y_ref_N);

        % solve the optimization problem of MPC at the current iteration
        ocp.solve();
        status = ocp.get('status');
        
        if status==0
            % no action if successfully solved
            % fprintf('\nsuccess!\n\n');
        else
            display('solution failed!');
        end
        
        % get solution
        x_traj = ocp.get('x');
        u_traj = ocp.get('u');
        pi_traj = ocp.get('pi');
        u_sim(:,ii) = ocp.get('u', 0);
        
        % shift trajectory for initialization
        x_traj_init = [x_traj(:,2:end), x_traj(:,end)];
        u_traj_init = [u_traj(:,2:end), u_traj(:,end)];
        pi_traj_init = [pi_traj(:,2:end), pi_traj(:,end)];
        
        % set current state of sim
        sim.set('x', x_sim(:,ii));
        % set current input in sim
        sim.set('u', u_sim(:,ii));
        
        % update theta_gradient with RMSE loss
        theta_gradient = theta_gradient + 2*([x_sim(1:2,ii)-desired_state(1:2,ii);zeros(1,1)])'*dx_dtheta;

        x_opt = ocp.get('x');
        u_opt = ocp.get('u');
    
        N = ocp_N;
        
        % Compute du_dxinit
        sens_u = zeros(nu, nx); % sens_u here is the du_dxinit, and we use built-in function in acados
        field = 'ex';
        stage = 0;
        % get sensitivities w.r.t. initial state value with index
        for index = 0:nx-1
            ocp.eval_param_sens(field, stage, index);
            temp = ocp.get('sens_u');
            sens_u(:,index+1) = temp(:,1);
        end
        du_dxinit = sens_u;
     
        % Solve for other analytical gradients (du/dtheta) by solving LMPC-Grads using quadprog
        du_dQR = get_dQR_quadprog(W,x_opt,u_opt,model,desired_state,ii,ocp_N,grad_f_x_fcn,grad_f_u_fcn);

        % use the current value of the state x and control input u to find
        % dx/dtheta following equation (5)
        dx_dtheta = (grad_f_x_fcn(x_sim(:,ii),u_sim(:,ii)) + grad_f_u_fcn(x_sim(:,ii),u_sim(:,ii))*du_dxinit)*dx_dtheta + grad_f_u_fcn(x_sim(:,ii),u_sim(:,ii))*du_dQR;

        % simulate the next state of the system based on the solution to
        % MPC optimziation
        sim.solve();
    
        % get new state
        x_sim(:,ii+1) = sim.get('xn');
    end
    u_hist = [u_hist; u_sim];
    % use RMSE as the loss
    RMSE_hist = [RMSE_hist sqrt(mean(sum((desired_state(1:2,1:1+T/model.Ts)-x_sim(1:2,:)).^2,1)))];
    loss_hist = [loss_hist sum((desired_state(1:2,1:1+T/model.Ts)-x_sim(1:2,:)).^2,'all')];
    theta_gradient_hist = [theta_gradient_hist;theta_gradient];
    fprintf('summed position error is %.3f\n',loss_hist(end));
    % update the cost coefficient
    W_new = W - learningRate * diag(theta_gradient);
    W_new = full(W_new);

    % sanity check
    coef_at_min = find(diag(W_new)<=minCoef);
    if ~isempty(coef_at_min)
        display('touching min coeff.')
        for jj = 1:length(coef_at_min)
            W_new(coef_at_min(jj),coef_at_min(jj)) = minCoef;
        end
    end

    coef_at_max = find(diag(W_new)>=maxCoef);
    if ~isempty(coef_at_max)
        display('touching max coeff.')
        for jj = 1:length(coef_at_max)
            W_new(coef_at_max(jj),coef_at_max(jj)) = maxCoef;
        end
    end
    
    W = W_new;
    % update the cost coefficient in the ocp
    ocp.set('cost_W',W_new);
    param_hist = [param_hist diag(W_new)];

    % compute tracking RMSE
    tracking_RMSE = sqrt(mean(sum((desired_state(1:2,1:1+T/model.Ts)-x_sim(1:2,:)).^2,1)));
    fprintf('RMSE is %.3f\n',tracking_RMSE);

    if difftune_itr >= total_itr
        break;
    else
        clf
    end

end

function [ocp,sim] = setup_acados(model,ocp_N,W)
    % handy arguments
    compile_interface = 'auto';
    codgen_model = 'true';
    % simulation
    sim_method = 'irk';
    sim_sens_forw = 'false';
    sim_num_stages = 4;
    sim_num_steps = 4;
    % ocp
    %ocp_nlp_solver = 'sqp';
    ocp_nlp_solver = 'sqp_rti';
    nlp_solver_max_iter = 100; % default
    ocp_qp_solver = 'partial_condensing_hpipm';
    %ocp_qp_solver = 'full_condensing_hpipm';
    ocp_qp_solver_cond_N = 5;
    ocp_sim_method = 'erk';
    %ocp_sim_method = 'irk';
    ocp_sim_method_num_stages = 2;
    ocp_sim_method_num_steps = 2;
    ocp_cost_type = 'linear_ls';
    %ocp_cost_type = 'nonlinear_ls';
    %ocp_cost_type = 'ext_cost';
    % Linear least square cost
    % linear least square cost: y^T * W * y, where y = Vx * x + Vu * u - y_ref
    nx = model.nx; % number of states
    nu = model.nu; % number of inputs
    ny = model.ny; % number of outputs in lagrange term
    ny_e = model.ny_e; % number of outputs in mayer term
    ng = model.ng;
    nh = model.nh;
    Vx = eye(ny, nx); % state-to-output matrix in lagrange term
    Vu = zeros(ny, nu);
    Vu(nx+1:end, :) = eye(nu); % input-to-output matrix in lagrange term
    Vx_e = Vx(1:ny_e,:); % state-to-output matrix in mayer term
    
    W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
    
    % constraints
    % based on the minimum and maximum angular velocity of each wheel, derive the constraint functions
    d = model.d;
    w_min = model.w_min;
    w_max = model.w_max;
    r = model.r;
    % follow acados definition: lg < Cx + Du < ug
    C = zeros(2,nx);
    D = [2 d; 2 -d];
    lg = 2*[w_min*r;w_min*r]; %lower bound
    ug = 2*[w_max*r;w_max*r]; %upper bound
    
    model.C = C;
    model.D = D;
    model.lg = lg;
    model.ug = ug;
    % acados ocp model
    ocp_model = acados_ocp_model();
    ocp_model.set('T', model.Ts*ocp_N); % set a much shorter horizon for MPC compared with total horizon
    
    % symbolics
    ocp_model.set('sym_x', model.sym_x);
    if isfield(model, 'sym_u')
        ocp_model.set('sym_u', model.sym_u);
    end
    if isfield(model, 'sym_xdot')
        ocp_model.set('sym_xdot', model.sym_xdot);
    end
    
    % cost
    ocp_model.set('cost_type', ocp_cost_type);
    ocp_model.set('cost_type_e', ocp_cost_type);
    if (strcmp(ocp_cost_type, 'linear_ls'))
        ocp_model.set('cost_Vu', Vu);
        ocp_model.set('cost_Vx', Vx);
        ocp_model.set('cost_Vx_e', Vx_e);
        ocp_model.set('cost_W', W);
        ocp_model.set('cost_W_e', W_e);
    elseif (strcmp(ocp_cost_type, 'nonlinear_ls'))
        ocp_model.set('cost_expr_y', model.expr_y);
        ocp_model.set('cost_expr_y_e', model.expr_y_e);
        ocp_model.set('cost_W', W);
        ocp_model.set('cost_W_e', W_e);
        ocp_model.set('cost_y_ref', yr);
        ocp_model.set('cost_y_ref_e', yr_e);
    else % if (strcmp(ocp_cost_type, 'ext_cost'))
        ocp_model.set('cost_expr_ext_cost', model.expr_ext_cost);
        ocp_model.set('cost_expr_ext_cost_e', model.expr_ext_cost_e);
    end
    
    % dynamics
    if (strcmp(ocp_sim_method, 'erk'))
        ocp_model.set('dyn_type', 'explicit');
        ocp_model.set('dyn_expr_f', model.expr_f_expl);
    else % irk
        ocp_model.set('dyn_type', 'implicit');
        ocp_model.set('dyn_expr_f', model.expr_f_impl);
    end
    
    % constraints
    x0 = zeros(nx, 1); %initial condition
    ocp_model.set('constr_x0', x0);
    if (ng>0)
	    ocp_model.set('constr_C', C);
	    ocp_model.set('constr_D', D);
	    ocp_model.set('constr_lg', lg);
	    ocp_model.set('constr_ug', ug);
    elseif (nh>0)
	    ocp_model.set('constr_expr_h', model.expr_h);
	    ocp_model.set('constr_lh', lh);
	    ocp_model.set('constr_uh', uh);
	    ocp_model.set('constr_expr_h_e', model.expr_h_e);
	    ocp_model.set('constr_lh_e', lh_e);
	    ocp_model.set('constr_uh_e', uh_e);
    else
	    ocp_model.set('constr_Jbx', Jbx);
	    ocp_model.set('constr_lbx', lbx);
	    ocp_model.set('constr_ubx', ubx);
	    ocp_model.set('constr_Jbu', Jbu);
	    ocp_model.set('constr_lbu', lbu);
	    ocp_model.set('constr_ubu', ubu);
    end
    
    % acados ocp opts
    ocp_opts = acados_ocp_opts();
    ocp_opts.set('compile_interface', compile_interface);
    ocp_opts.set('codgen_model', codgen_model);
    ocp_opts.set('param_scheme_N', ocp_N);
    ocp_opts.set('nlp_solver', ocp_nlp_solver);
    ocp_opts.set('qp_solver', ocp_qp_solver);
    ocp_opts.set('nlp_solver_max_iter', nlp_solver_max_iter);
    if (strcmp(ocp_qp_solver, 'partial_condensing_hpipm'))
        ocp_opts.set('qp_solver_cond_N', ocp_qp_solver_cond_N);
    end
    ocp_opts.set('sim_method', ocp_sim_method);
    ocp_opts.set('sim_method_num_stages', ocp_sim_method_num_stages);
    ocp_opts.set('sim_method_num_steps', ocp_sim_method_num_steps);
    ocp_opts.set('regularize_method', 'no_regularize');
    
    % acados ocp (compiling mex files)
    % create ocp
    ocp = acados_ocp(ocp_model, ocp_opts);
    
    % acados sim model
    sim_model = acados_sim_model();
    % symbolics
    sim_model.set('sym_x', model.sym_x);
    if isfield(model, 'sym_u')
        sim_model.set('sym_u', model.sym_u);
    end
    if isfield(model, 'sym_xdot')
        sim_model.set('sym_xdot', model.sym_xdot);
    end
    
    % model
    sim_model.set('T', model.Ts);
    if (strcmp(sim_method, 'erk'))
        sim_model.set('dyn_type', 'explicit');
        sim_model.set('dyn_expr_f', model.expr_f_expl);
    else % irk
        sim_model.set('dyn_type', 'implicit');
        sim_model.set('dyn_expr_f', model.expr_f_impl);
    end
    
    % acados sim opts
    sim_opts = acados_sim_opts();
    sim_opts.set('compile_interface', compile_interface);
    sim_opts.set('codgen_model', codgen_model);
    sim_opts.set('num_stages', sim_num_stages);
    sim_opts.set('num_steps', sim_num_steps);
    sim_opts.set('method', sim_method);
    sim_opts.set('sens_forw', sim_sens_forw);
    
    % acados sim (compiling mex)
    % create sim
    sim = acados_sim(sim_model, sim_opts);
    % sim.C_sim
    % sim.C_sim_ext_fun
end

function du_dQR = get_dQR_quadprog(W,x_opt,u_opt,model,desired_state,ii,ocp_N,grad_f_x_fcn,grad_f_u_fcn)
    % use quadprog to compute the gradient
    du_dQR = [];
    N = ocp_N;
    nx = model.nx;
    nu = model.nu;
    Q = W(1:nx,1:nx);
    R = W(nx+1:end,nx+1:end);

    % for the data loading below, we use the notation with quadprog
    H = [];
    Aeq = [];
    Beq = [];
    for kk = 1:N
        % form the H matrix
        H = blkdiag(H,W);
        
        % form Aeq and Beq
        if kk == 1
            Aeq = [eye(nx) zeros(nx,N*(nx+nu))];
        else
            A = full(grad_f_x_fcn(x_opt(:,kk-1),u_opt(:,kk-1)));
            B = full(grad_f_u_fcn(x_opt(:,kk-1),u_opt(:,kk-1)));
            Aeq = [Aeq;
                 zeros(nx,(nx+nu)*(kk-2)) [-A -B eye(nx)] zeros(nx,(nx+nu)*(N-kk+1))];
        end
        % Beq is a zero vecotr
    end
    A = full(grad_f_x_fcn(x_opt(:,N),u_opt(:,N)));
    B = full(grad_f_u_fcn(x_opt(:,N),u_opt(:,N)));
    Aeq = [Aeq;zeros(nx,(nx+nu)*(N-1)) [-A -B eye(nx)] ];
    Beq = zeros(size(Aeq,1),1);
    H = blkdiag(H,Q);

    % only f needs to be changed everytime
    identity = eye(nx+nu);
    % jj represents the elements in W
    for jj = 1:nx+nu
        unitVec = identity(:,jj);
        f = [];
        for kk = 1:N
            % form f vector
            tau_opt = [x_opt(:,kk);u_opt(:,kk)];
            tau_ref = [desired_state(:,ii+kk-1);
                zeros(nu,1)];
            
            f = [f;[tau_opt(jj)-tau_ref(jj)]*unitVec]; 
        end
        tau_opt = [x_opt(:,N+1)];
        tau_ref = [desired_state(:,ii+N)];
        % handle the tail cases
        if jj<=nx
            f = [f;(tau_opt(jj)-tau_ref(jj))*unitVec(1:nx)]; 
        else
            f = [f;0*unitVec(1:nx)]; 
        end
        grad = quadprog(H,f,[],[],Aeq,Beq,[],[],[]);
        du_dQR = [du_dQR grad(nx+1:nx+nu)]; 
    end
end


