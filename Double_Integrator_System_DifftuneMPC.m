% This script implements Difftune-MPC on a double integrator system.
% The original MPC problem is solved using acaods and the analytical
% gradients are obtained by solving extra MPC problems (LMPC-Grad) using
% acados.
% Please install acados in matlab before running this example.
% We define constraints, cost functions and everything in this example 
% following instructions from acados.
% For more info, check problem_formulation_ocp_mex.pdf from acados

% Ran Tao, Sheng Cheng
% University of Illinois Urbana-Champaign

clear all
clc
close all

%% setup problem
model = Double_Integrator_System;

% dims
T = 10.0; % total horizon time for the simulation
nx = model.nx; % number of states
nu = model.nu; % number of inputs
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
nbx = nx; % number of state bounds
nbu = nu; % number of input bounds
u_max = model.u_max;
x_max = model.x_max;
v_max = model.v_max;

%% Set up acados
ocp_N = 20; % prediction horizon for discrete time system in MPC
% Initial weight matrix in lagrange term
W = diag([1*ones(nx,1);... % pos, vel
            1*ones(nu,1) ]); % control
Jbx = eye(nbx, nx); 
lbx = -[x_max;v_max]; %lower state bound
ubx =  [x_max;v_max]; %upper state bound
Jbu = eye(nbu, nu); 
lbu = -u_max; %lower control input bound
ubu =  u_max; %upper control input bound
[ocp,sim] = setup_acados(model,ocp_N,W,Jbx,lbx,ubx,Jbu,lbu,ubu);
ocp_grad=ocp; %ocp is for solving original MPC problem and ocp_grad is for solving LMPC-Grad.

%% add desired trajectory for tracking
pos_des = @(t) 1*t + 1*(1 - cos(t));
vel_des = @(t) 1*ones(size(t)) + 1*sin(t);

%% Difftune
% Parameters for DiffTune
learningRate = 0.01; 
minCoef = 0.01; % the minimial diagonal elements of Q and R cannot be smaller than this value
maxCoef = 100; % the minimial diagonal elements of Q and R cannot be bigger than this value
loss_hist = [];
param_hist = [diag(W)];
RMSE_hist = [];

% initialize the gradient
theta_gradient = zeros(1,nx+nu);
theta_gradient_hist = theta_gradient;
total_itr = 20;
difftune_itr = 0;

%% DiffTune main loop
while(1)
    difftune_itr = difftune_itr + 1;

    %% closed loop simulation
    n_sim = T/model.Ts;
    x_sim = zeros(nx, n_sim+1);
    x_sim(:,1) = zeros(nx,1);
    u_sim = zeros(nu, n_sim);
    
    % set the initial guess of x, u, and pi (costate)
    x_traj_init = zeros(nx, ocp_N+1);
    u_traj_init = zeros(nu, ocp_N);
    pi_traj_init = zeros(nx, ocp_N);

    % desired state
    desired_traj = pos_des(model.Ts * [0:n_sim]);
    desired_state = [desired_traj;
        vel_des(model.Ts * [0:n_sim])];

    % initialize the sensitivity
    dx_dtheta = zeros(nx,nx+nu);
    theta_gradient = zeros(1,nx+nu);

    % for each iteration in the following for-loop, we solve the MPC 
    % problem using acados, and compute the Jacobians dx/dtheta and 
    % du/dtheta in equation (5).
    % After the for-loop, we have run simulation over entire time hoziron 
    % T, and then use the calculated dx/dtheta and du/dtheta to update 
    % theta, which contains Q and R matrices that are compactly coded in 
    % the variable W in this code.

    for ii=1:n_sim
        % reset the initial states and references
        ocp.set('constr_x0', x_sim(:,ii));
        ocp.set('init_x', x_traj_init);
        ocp.set('init_u', u_traj_init);
        ocp.set('init_pi', pi_traj_init);

        % set the y_ref and y_ref_e for tracking a given signal
        time_varying_y_ref = [pos_des((ii+[0:ocp_N-1])*model.Ts);
            vel_des((ii+[0:ocp_N-1])*model.Ts);
            zeros(1,ocp_N)];
        for jj = 1:ocp_N
            ocp.set('cost_y_ref', time_varying_y_ref(:,jj), jj-1);
        end
        time_varying_y_ref_N = [pos_des((ii+ocp_N)*model.Ts);
            vel_des((ii+ocp_N)*model.Ts)];
        ocp.set('cost_y_ref_e', time_varying_y_ref_N);

        % solve the optimization problem of MPC at the current step
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

        % set initial state of sim
        sim.set('x', x_sim(:,ii));
        % set input in sim
        sim.set('u', u_sim(:,ii));

        % update theta_gradient, which is corresponding to dL/dtheta in the
        % paper
        theta_gradient = theta_gradient + 2*([x_sim(1,ii)-desired_state(1,ii);0])'*dx_dtheta;

        % get du/dQR (diag), with analytical grad
        x_opt = ocp.get('x');
        u_opt = ocp.get('u');

        N = ocp_N;

        if abs(u_opt(1)-ubu)<=0.01 || abs(u_opt(1)-lbu)<=0.01 
            % if the control constraints are active, then analytical 
            % gradients will be zero
            du_dQR = [0 0 0];
            du_dxinit = [0 0];
        else
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
            % Ran, can you wrap the lines of code below to get_dQR_acados?
            % similar to how you wrap the function get_dQR_quadprog in the
            % other script.
            % Compute du_dQ and du_dR by solving LMPC-Grad using ocp-grad
            ocp_grad.set('constr_x0', zeros(nx,1));

            % reformulate the cost function to solve for sensitivity
            RHS = zeros(N*(nx+nu),(nx+nu)^2);
            for kk = 1:N
                tau_opt = [x_opt(:,kk);u_opt(:,kk)];
                tau_ref = [pos_des((ii+kk-1)*model.Ts);
                    vel_des((ii+kk-1)*model.Ts);
                    zeros(nu,1)];
                RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,:) = 0.5* (kron((-tau_opt+tau_ref)',eye(nx+nu)) + kron(eye(nx+nu),(-tau_opt+tau_ref)'));
            end

            % solve for the du/dQ11
            for kk = 1:N
                ocp_grad.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,1),kk-1); % set the 1st column of RHS, which correspond to Q11
            end
            ocp_grad.set('cost_y_ref_e', zeros(nx,1));
            ocp_grad.solve;
            du_dQ11 = ocp_grad.get('u',0);

            % solve for the du/dQ22
            for kk = 1:N
                ocp_grad.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,5),kk-1); % set the 5th column of RHS, which correspond to Q22
            end
            ocp_grad.set('cost_y_ref_e', zeros(nx,1));
            ocp_grad.solve;
            du_dQ22 = ocp_grad.get('u',0);

            % solve for the du/dR
            for kk = 1:N
                ocp_grad.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,9),kk-1); % set the 9th column of RHS, which correspond to R
            end
            ocp_grad.set('cost_y_ref_e', zeros(nx,1));
            ocp_grad.solve;
            du_dR = ocp_grad.get('u',0);

            du_dQR = [du_dQ11 du_dQ22 du_dR];
        end

        % sensitivity propagation: use the current value of the state x and 
        % control input u to find dx/dtheta following equation (5)
        dx_dtheta = (model.A + model.B*du_dxinit)*dx_dtheta + model.B*du_dQR;

        % simulate state
        sim.solve();

        % get new state
        x_sim(:,ii+1) = sim.get('xn');
    end
    % save the history of RMSE
    RMSE_hist = [RMSE_hist sqrt(mean(sum((desired_state(1,:)-x_sim(1,:)).^2,1)))];
    % save the history of the loss
    loss_hist = [loss_hist sum((desired_state(1,:)-x_sim(1,:)).^2,'all')];
    % save the history of theta_gradient
    theta_gradient_hist = [theta_gradient_hist;theta_gradient];

    fprintf('summed position error is %.3f\n',loss_hist(end));

    % update the cost coefficient
    W_new = W - learningRate * diag(theta_gradient);

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
    ocp_grad.set('cost_W',W_new);
    param_hist = [param_hist diag(W_new)];

    % compute tracking RMSE
    tracking_RMSE = sqrt(mean((desired_traj-x_sim(1,:)).^2));
    fprintf('RMSE is %.3f\n',tracking_RMSE);

    if difftune_itr == 1
        x_sim_init = x_sim;
    end

    if difftune_itr == total_itr
        x_sim_final = x_sim;
    end

    if difftune_itr >= total_itr
        break;
    end
end

% tracking performance
figure;
subplot(2,1,1);
plot((0:n_sim)*model.Ts, x_sim_init(1,:),'r-.','linewidth',2,'DisplayName','w. initial parameters');
hold on;
plot((0:n_sim)*model.Ts, x_sim_final(1,:),'b-','linewidth',2,'DisplayName','w. learned parameters');
plot((0:n_sim)*model.Ts,desired_state(1,1:n_sim+1),'-','color',[0 0 1 0.4],'linewidth',2,'DisplayName','desired position');
xlabel('Time [s]');
ylabel('x [m]');
legend('Location','best');
title('position tracking comparison between initial and learned parameters');

subplot(2,1,2);
plot((0:n_sim)*model.Ts, x_sim_init(2,:),'r-.','linewidth',2,'DisplayName','w. initial parameters');
hold on;
plot((0:n_sim)*model.Ts, x_sim_final(2,:),'b-','linewidth',2,'DisplayName','w. learned parameters');
plot((0:n_sim)*model.Ts,desired_state(2,1:n_sim+1),'-','color',[0 0 1 0.4],'linewidth',2,'DisplayName','desired position');
xlabel('Time [s]');
ylabel('y [m]');
legend('Location','northeast');

% RMSE reduction
figure;
plot(RMSE_hist,'LineWidth',2);
xlabel('iterations');
ylabel('RMSE [m]');


function [ocp,sim] = setup_acados(model,ocp_N,W,Jbx,lbx,ubx,Jbu,lbu,ubu)
    % handy arguments
    compile_interface = 'auto';
    codgen_model = 'true';
    % simulation
    sim_method = 'irk';
    sim_sens_forw = 'false';
    sim_num_stages = 4;
    sim_num_steps = 4;
    % ocp
    ocp_nlp_solver = 'sqp';
    %ocp_nlp_solver = 'sqp_rti';
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
    ny = nu+nx; % number of outputs in lagrange term
    ny_e = nx; % number of outputs in mayer term
    Vx = eye(ny, nx); % state-to-output matrix in lagrange term
    Vu = zeros(ny, nu);
    Vu(nx+1:end, :) = eye(nu); % input-to-output matrix in lagrange term
    Vx_e = Vx(1:ny_e,:); % state-to-output matrix in mayer term
    W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term
    
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
    x0 = zeros(nx, 1); %initial condition, redefined later
    ocp_model.set('constr_x0', x0);
    ocp_model.set('constr_Jbx', Jbx);
    ocp_model.set('constr_lbx', lbx);
    ocp_model.set('constr_ubx', ubx);
    ocp_model.set('constr_Jbu', Jbu);
    ocp_model.set('constr_lbu', lbu);
    ocp_model.set('constr_ubu', ubu);
    
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