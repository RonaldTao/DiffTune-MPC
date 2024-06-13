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
% This script implements difftune for MPC
% The analytical gradients are obtained by solving extra MPC problems
% (rather than obtained by KKT factorization)

% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
    ! source env.sh
    % 	error('env.sh has not been sourced! Before executing this example, run: source env.sh');
end

% for DiffTune
learningRate = 0.01; % 0.5 works for ramp trajectory
minCoef = 0.01; % the minimial diagonal elements of Q and R cannot be smaller than this value
maxCoef = 1000; % the minimial diagonal elements of Q and R cannot be bigger than this value
%% handy arguments
compile_interface = 'auto';
codgen_model = 'true';
% simulation
sim_method = 'irk';
sim_sens_forw = 'false';
sim_num_stages = 4;
sim_num_steps = 4;
% ocp
ocp_N = 10; % N is the prediction horizon
%ocp_nlp_solver = 'sqp';
ocp_nlp_solver = 'sqp_rti';
nlp_solver_max_iter = 100; % default
ocp_qp_solver = 'partial_condensing_hpipm';
%ocp_qp_solver = 'full_condensing_hpipm';
ocp_qp_solver_cond_N = 5;
ocp_sim_method = 'erk'; % sheng: change to erk to make simulation faster
%ocp_sim_method = 'irk';
ocp_sim_method_num_stages = 2;
ocp_sim_method_num_steps = 2;
ocp_cost_type = 'linear_ls';
%ocp_cost_type = 'nonlinear_ls';
%ocp_cost_type = 'ext_cost';

%% setup problem
% linear mass-spring system
model = Unicycle_model;
%f=get_jacobian(model);
% set the reference control to be always 0
uref = 0;


% dims
T = 10.0; % horizon length time in seconds
nx = model.nx; % number of states
nu = model.nu; % number of inputs
ny = nu+nx; % number of outputs in lagrange term
ny_e = nx; % number of outputs in mayer term
nbx = nx; % number of state bounds
nbu = nu; % number of input bounds

% constraint formulation
% bounds on x and u
nbx = nx;
nbu = nu;
ng = 1;
nh = 0;

%% we use linear cost here
% see the official guide here: file:///home/sheng/Downloads/problem_formulation_ocp_mex.pdf
% see the example from
% /home/sheng/acados/examples/acados_matlab_octave/pendulum_dae/example_closed_loop.m
% linear least square cost: y^T * W * y, where y = Vx * x + Vu * u - y_ref
Vx = eye(ny, nx); % state-to-output matrix in lagrange term

Vu = zeros(ny, nu);
Vu(nx+1:end, :) = eye(nu); % input-to-output matrix in lagrange term

Vx_e = Vx(1:ny_e,:); % state-to-output matrix in mayer term

% weight matrix in lagrange term
W = diag([1*ones(nx,1);... % pos, vel
            1*ones(nu,1) ]); % control
W_e = W(nu+1:nu+nx, nu+1:nu+nx); % weight matrix in mayer term

% Sheng: set the yr and yr_e in a tracking control fashion
% yr = [xtarget; uref]; % output reference in lagrange term
% yr_e = xtarget; % output reference in mayer term

% constraints
x0 = zeros(nx, 1);
% Sheng: The prefix J (e.g. Jbx) represents a matrix made of rows from an
%        identity matrix. This is a way to encode which of your states or
%        inputs are bounded.
% use dummy bound for solver to proceed without error
Jbx = eye(nbx, nx); % for ii=1:nbx Jbx(ii,ii)=1.0; end
lbx = -100*ones(nbx, 1);
ubx =  100*ones(nbx, 1);
Jbu = eye(nbu, nu); % for ii=1:nbu Jbu(ii,ii)=1.0; end
lbu = -100*ones(nu, 1);
ubu =  100*ones(nu, 1);
d=model.d;
w_min = model.w_min;
w_max = model.w_max;
r = model.r;
% set bounds on the total control actions such that
% u1+u2<= 1
% u1+u2>=-1

% use C x(t) + Du(t) − Jsg su,g ≤ bar{g} in acados
%     \lowerbar{g} ≤ C x(t) + D u(t) + Jsg sl,g(t),
C = zeros(2,nx);

D = zeros(2,nu);
D = [2 d; 2 -d];

lg = 2*[w_min*r;w_min*r];
ug = 2*[w_max*r;w_max*r];

model.C = C;
model.D = D;
model.lg = lg;
model.ug = ug;
%% acados ocp model
ocp_model = acados_ocp_model();
ocp_model.set('T', model.Ts*ocp_N); % set a much shorter horizon for MPC
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
    % 	ocp_model.set('cost_y_ref', yr);
    % 	ocp_model.set('cost_y_ref_e', yr_e);
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
ocp_model.set('constr_x0', x0);
if (ng>0)
	ocp_model.set('constr_C', C);
	ocp_model.set('constr_D', D);
	ocp_model.set('constr_lg', lg);
	ocp_model.set('constr_ug', ug);
    % we do not need the following path constraints on the temrinal state
    % for now
% 	ocp_model.set('constr_C_e', C_e);
% 	ocp_model.set('constr_lg_e', lg_e);
% 	ocp_model.set('constr_ug_e', ug_e);
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

%% acados ocp opts
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

%% acados ocp (compiling mex files)
% create ocp
ocp = acados_ocp(ocp_model, ocp_opts);
% ocp.C_ocp
% ocp.C_ocp_ext_fun

%% acados sim model
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
% Note: "T" here stands for the time of simulation for one iteration/sample period
% The total time of simulation is determined by n_sim, which defines the
% number of samples/steps of a discrete-time simulation
if (strcmp(sim_method, 'erk'))
    sim_model.set('dyn_type', 'explicit');
    sim_model.set('dyn_expr_f', model.expr_f_expl);
else % irk
    sim_model.set('dyn_type', 'implicit');
    sim_model.set('dyn_expr_f', model.expr_f_impl);
end

%% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('compile_interface', compile_interface);
sim_opts.set('codgen_model', codgen_model);
sim_opts.set('num_stages', sim_num_stages);
sim_opts.set('num_steps', sim_num_steps);
sim_opts.set('method', sim_method);
sim_opts.set('sens_forw', sim_sens_forw);

%% acados sim (compiling mex)
% create sim
sim = acados_sim(sim_model, sim_opts);
% sim.C_sim
% sim.C_sim_ext_fun

%% add desired trajectory for tracking
%% add desired trajectory for tracking
pos_1_des = @(t) 1-cos(t/2);
pos_2_des = @(t) 0.5*t;
pos_3_des = @(t) atan2(0.5,0.5*sin(t/2));
%% Difftune
loss_hist = [];
param_hist = [diag(W)];
RMSE_hist = [];
theta_hist = [];
u_hist = [];
% initialize the gradientF
theta_gradient = zeros(1,nx+nu); 
theta_gradient_hist = theta_gradient;
total_itr=20;
difftune_itr = 0;

import casadi.*
grad_f_x=jacobian(model.expr_f_discrete,model.sym_x);
grad_f_x_fcn = Function('grad_f_x_fcn',{model.sym_x,model.sym_u},{grad_f_x});

grad_f_u=jacobian(model.expr_f_discrete,model.sym_u);
grad_f_u_fcn = Function('grad_f_u_fcn',{model.sym_x,model.sym_u},{grad_f_u});






% desired state
%% loop
while(1)
    theta_gradient = zeros(1,nx+nu);
    difftune_itr = difftune_itr + 1;
    %% closed loop simulation
    n_sim = T/model.Ts;
    x_sim = zeros(nx, n_sim+1);
    x_sim(:,1) = [0;0;pi/2];
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
    % for each iteration in the following loop, we solve the MPC problem, and achieve the
    % dx/dtheta and du/dtheta
    % After the loop, we have run simulation over entire time hoziron T, and then use
    % the calculated dx/dtheta and du/dtheta to update theta, which is Q
    % and R in W
    
    timing_ocp = [];
    timing_grad = [];
    timing_grad2 = [];
    ocp_grad = ocp;
    tic;
    error_percent = [];
    for ii=1:n_sim
        ocp.model_struct.dyn_type = 'explicit';
        ocp.model_struct.dyn_expr_f = model.expr_f_expl;
        % reset constraints back to normal
        ocp.set('constr_C',C);
        ocp.set('constr_D',D);
        ocp.set('constr_ug',ug);
        ocp.set('constr_lg',lg);

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

        
        ocp_tic = tic;

        % solve the optimization problem of MPC at the current iteration

        ocp.solve();
        ii;
        timing_ocp = [timing_ocp toc(ocp_tic)];
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

        u_traj;
        t = D*u_traj;
        for j = 1:ocp_N*2
            if t(j)>ug
                1;
                t;
            elseif t(j)<lg
                2;
                t;
            end
        end
        ocp.model_struct.constr_lg;
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
        
        % update theta_gradient with RMSE loss (we will do the scaling until the end)
    %         theta_gradient = theta_gradient + (x_sim(:,ii)-desired_state(:,ii))'*dx_dtheta;
            
        % update theta_gradient with sumeed loss (no scaling necessary)
        % use total state only
    %         theta_gradient = theta_gradient + 2*([x_sim(:,ii)-desired_state(:,ii)])'*dx_dtheta;
        % use position error only equation (4) of one step
        theta_gradient = theta_gradient + 2*([x_sim(1:2,ii)-desired_state(1:2,ii);zeros(1,1)])'*dx_dtheta;
        % if ii<10
        %     ii
        %     ocp.model_struct.dyn_expr_f
        %     dx_dtheta
        %     current_pos_error=([x_sim(1:3+4,ii)-desired_state(1:3+4,ii);zeros(3+3,1)])
        %     %theta_gradient
        % 
        % end
        % get du/dQR (diag), with analytical grad
        x_opt = ocp.get('x');
        u_opt = ocp.get('u');
    
        N = ocp_N;
        
        % comparison with acados's sensitivity to initial state
        sens_u = zeros(nu, nx); % storage of du_dxinit of dim nu x nx
        field = 'ex'; % equality constraint on states
        stage = 0;
        % get sensitivities w.r.t. initial state value with index
        for index = 0:nx-1
            ocp.eval_param_sens(field, stage, index);
            temp = ocp.get('sens_u');
            sens_u(:,index+1) = temp(:,1);
        end
    %         disp('solution sensitivity dU_dx0')
    %         disp(sens_u)
    
        du_dxinit = sens_u; % Sheng: start here

        % set the active constrains only for the inequality constraints
        C = model.C;
        D = model.D;
        lg = model.lg;
        ug = model.ug;

        % enforce the constraints for gradient computation in an alternative way
        % by reforming C and D

%         for kk = 1:N
%             idx_act = find(abs(C*x_opt(:,kk)+D*u_opt(:,kk) - ug)<1e-5 | abs(C*x_opt(:,kk)+D*u_opt(:,kk) - lg)<1e-5);
%             new_C = zeros(size(C));
%             new_D = zeros(size(D));
%             
%             if ~isempty(idx_act)
%                 % if saturation happens on upper/lower bound, then only force the
%                 % corresponding rows of C and D to be same as the original problem
%                 new_C(idx_act,:) = C(idx_act,:);
%                 new_D(idx_act,:) = D(idx_act,:);
%             else
%                 % no active constraints, set C and D to zero
%                 % no action needed
%             end
% 
%             ocp.set('constr_C', new_C, kk-1);
%             ocp.set('constr_D', new_D, kk-1);
%         end
%         % need to enforece the bounds to be zeros
%         act_const = zeros(length(ug),1);
%         ocp.set('constr_ug', act_const);
%         ocp.set('constr_lg', -act_const);

%         active_index = [];
%         for kk = 1:N
%             idx_act = find(abs([C D;-C -D]*[x_opt(:,kk);u_opt(:,kk)] - [ug;-lg])<1e-5);
%             if ~isempty(idx_act)
%                 % if saturation happens on upper/lower bound, then only force the
%                 % corresponding rows of C and D to be same as the original problem
%                 active_index = [active_index idx_act];
%             else
%                 active_index = [active_index 0];
%                 % no active constraints, set C and D to zero
%                 % no action needed
%             end
%         end

        grad_f_x2d=jacobian(model.expr_f_discrete,model.sym_x);
        grad_f_x_fcn2d = Function('grad_f_x_fcn2d',{model.sym_x,model.sym_u},{grad_f_x2d});
    
        grad_f_u2d=jacobian(model.expr_f_discrete,model.sym_u);
        grad_f_u_fcn2d = Function('grad_f_u_fcn2d',{model.sym_x,model.sym_u},{grad_f_u2d});
        
        grad_f_x2c=jacobian(model.expr_f_expl,model.sym_x);
        grad_f_x_fcn2c = Function('grad_f_x_fcn2c',{model.sym_x,model.sym_u},{grad_f_x2c});
    
        grad_f_u2c=jacobian(model.expr_f_expl,model.sym_u);
        grad_f_u_fcn2c = Function('grad_f_u_fcn2d',{model.sym_x,model.sym_u},{grad_f_u2c});
        % For ocp, the only different componenet will be the new
        % linearlized dynamics model
        %du_dQR=get_dQR(x_opt,u_opt,model,desired_state,ii,ocp,ocp_N);
        grad_tic = tic;
        [du_dQR,timing_an] = get_dQR_quadprog(W,x_opt,u_opt,ocp_grad,model,desired_state,ii,ocp_N,grad_f_x_fcn2d,grad_f_u_fcn2d);
        timing_grad = [timing_grad toc(grad_tic)];
        
        grad_tic2 = tic;
        [du_dQR2,timing_nu] = get_dQR_numerical(W,nx,nu,ocp_grad,x_sim,x_traj_init,u_traj_init,pi_traj_init,time_varying_y_ref,desired_state,ii,ocp_N);
        timing_grad2 = [timing_grad2 toc(grad_tic2)];
        error_percent = [error_percent abs((du_dQR-du_dQR2)./du_dQR)];
        %find the gradient function of the nonlinear dynamics to find the
        %linearzied system using discrete time model
        % use the current value of the state x and control input u to find
        % dx/dtheta equation (5a)
        dx_dtheta = (grad_f_x_fcn(x_sim(:,ii),u_sim(:,ii)) + grad_f_u_fcn(x_sim(:,ii),u_sim(:,ii))*du_dxinit)*dx_dtheta + grad_f_u_fcn(x_sim(:,ii),u_sim(:,ii))*du_dQR;
        % if ii<10
        %     du_dQR
        %     current_x=x_sim(:,ii)
        %     current_u=u_sim(:,ii)
        %     current_dfx=full(grad_f_x_fcn(x_sim(:,ii),u_sim(:,ii)))
        %     du_dxinit
        %     current_dfu=full(grad_f_u_fcn(x_sim(:,ii),u_sim(:,ii)))
        %     ocp.model_struct.dyn_expr_f
        % end
        % simulate the next state of the system based on the solution to
        % MPC optimziation
        sim.solve();
    
        % get new state
        x_sim(:,ii+1) = sim.get('xn');
    end

    avg_time_solve = toc/n_sim

    if difftune_itr == 1
        save x_i.mat x_sim;
        save u_i.mat u_sim;
    end

    if difftune_itr == total_itr
        save x_f.mat x_sim;
        save u_f.mat u_sim;
        save ref.mat desired_state;
    end
    u_hist = [u_hist; u_sim];
    % use RMSE as the loss
    RMSE_hist = [RMSE_hist sqrt(mean(sum((desired_state(1:2,1:1+T/model.Ts)-x_sim(1:2,:)).^2,1)))];

    loss_hist = [loss_hist sum((desired_state(1:2,1:1+T/model.Ts)-x_sim(1:2,:)).^2,'all')];
    theta_gradient_hist = [theta_gradient_hist;theta_gradient];
    fprintf('summed position error is %.3f\n',loss_hist(end));
    % update the cost coefficient equation (3)
    W_new = W - learningRate * diag(theta_gradient);
    W_new = full(W_new);
    % % keep R22 as unity for reference of scaling
    % W_new(6,6) = W(6,6);

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
    %ocp_grad.set('cost_W',W_new);
    param_hist = [param_hist diag(W_new)];

    % compute tracking RMSE
    tracking_RMSE = sqrt(mean(sum((desired_state(1:2,1:1+T/model.Ts)-x_sim(1:2,:)).^2,1)));
    fprintf('RMSE is %.3f\n',tracking_RMSE);

    if difftune_itr >= total_itr
        break;
    else
        clf
    end

    if difftune_itr == 1
        figure;
    end

    set(gcf,'Position',[551 194 1019 582]);
    subplot(3, 3, [1:2])
    plot((0:n_sim)*model.Ts, x_sim(1,:),'b-','linewidth',2,'DisplayName','closed-loop-x');
    hold on;
    plot((0:n_sim)*model.Ts,desired_state(1,1:n_sim+1),'r-','linewidth',2,'DisplayName','desired-x');
    plot((0:n_sim)*model.Ts, x_sim(2,:),'b--','linewidth',2,'DisplayName','closed-loop-y');
    plot((0:n_sim)*model.Ts,desired_state(2,1:n_sim+1),'r--','linewidth',2,'DisplayName','desired-y');
%     plot((0:n_sim)*model.Ts, x_sim(3,:),'b:','linewidth',2,'DisplayName','closed-loop-z');
%     plot((0:n_sim)*model.Ts,desired_state(3,1:n_sim+1),'r:','linewidth',2,'DisplayName','desired-z');
    
    
    title(['closed loop simulation, trial=' num2str(difftune_itr)]);
    xlabel('time [s]');
    ylabel('pos [m]');
    legend
    subplot(3, 3, [4:5]);
    plot((0:n_sim)*model.Ts, x_sim(3,:)/pi*180,'b-','linewidth',2,'DisplayName','closed-loop-x');
    hold on;
    plot((0:n_sim)*model.Ts,desired_state(3,1:n_sim+1)/pi*180,'r-','linewidth',2,'DisplayName','desired-theta');
%     plot((0:n_sim)*model.Ts, x_sim(9,:),'b--','linewidth',2,'DisplayName','closed-loop-y');
%     plot((0:n_sim)*model.Ts,desired_state(9,1:n_sim+1),'r--','linewidth',2,'DisplayName','desired-y');
%     plot((0:n_sim)*model.Ts, x_sim(10,:),'b:','linewidth',2,'DisplayName','closed-loop-y');
%     plot((0:n_sim)*model.Ts,desired_state(10,1:n_sim+1),'r:','linewidth',2,'DisplayName','desired-y');

    xlabel('time [s]');
    ylabel('theta [degree]');

    subplot(3, 3, [7:8])
    plot((1:n_sim)*model.Ts, u_sim,'linewidth',2);
    ylim([-0.8 0.8]);
    ylabel('control u')
    xlabel('time [s]');

    subplot(3,3,[3;6;9]);
    plot(1:difftune_itr,RMSE_hist(1:difftune_itr),'linewidth',2);
    hold on;
    stem(difftune_itr,RMSE_hist(difftune_itr),'Color','b');
    xlabel('trial');
    ylabel('RMSE [m]');
    xlim([1 total_itr]);
    if difftune_itr<total_itr
        drawnow;
    end
    saveas(gcf,['unicycle_difftune' num2str(difftune_itr) '.png']);
      
end
save RMSE.mat RMSE_hist;
save u_hist.mat u_hist;
% error = error_percent(:);
% error = sort(error);
% mean(error);
% median(error);
%%
%% plot result
figure;
plot(loss_hist);
xlabel('iteration');
ylabel('loss');
title('loss reduction, analy grad with no KKT')

function [du_dQR,timing_an] = get_dQR_quadprog(W,x_opt,u_opt,ocp_grad,model,desired_state,ii,ocp_N,grad_f_x_fcn2d,grad_f_u_fcn2d)
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
    time = [];
    for kk = 1:N
        % form the H matrix
        H = blkdiag(H,W);
        
        % form Aeq and Beq
        if kk == 1
            Aeq = [eye(nx) zeros(nx,N*(nx+nu))];
        else
            A = full(grad_f_x_fcn2d(x_opt(:,kk-1),u_opt(:,kk-1)));
            B = full(grad_f_u_fcn2d(x_opt(:,kk-1),u_opt(:,kk-1)));
            Aeq = [Aeq;
                 zeros(nx,(nx+nu)*(kk-2)) [-A -B eye(nx)] zeros(nx,(nx+nu)*(N-kk+1))];
        end
        % Beq is a zero vecotr
    end
    A = full(grad_f_x_fcn2d(x_opt(:,N),u_opt(:,N)));
    B = full(grad_f_u_fcn2d(x_opt(:,N),u_opt(:,N)));
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
        tim = tic;
        grad = quadprog(H,f,[],[],Aeq,Beq,[],[],[]);
        time = [time toc(tim)];
        du_dQR = [du_dQR grad(nx+1:nx+nu)]; 
    end
    timing_an = sum(time);
end



function du_dQR=get_dQR(x_opt,u_opt,model,desired_state,ii,ocp,ocp_N) 
    W = ocp.model_struct.cost_W;
    
    Ac = grad_f_x_fcn2c(x_opt(:,1),u_opt(:,1));
    Bc = grad_f_u_fcn2c(x_opt(:,1),u_opt(:,1));

    expr_f_expl = Ac*model.sym_x + Bc*model.sym_u;
    expr_f_impl = expr_f_expl - model.sym_xdot;

    Ad = grad_f_x_fcn2d(x_opt(:,1),u_opt(:,1));
    Bd = grad_f_u_fcn2d(x_opt(:,1),u_opt(:,1));

    expr_phi = Ad*model.sym_x + Bd*model.sym_u;
    
    ocp.model_struct.dyn_type = 'discrete';
    ocp.model_struct.dyn_expr_f = expr_phi;

    % ocp.model_struct.dyn_type = 'explicit';
    % ocp.model_struct.dyn_expr_f = expr_f_expl;
    N = ocp_N;
    nx = model.nx;
    nu = model.nu;
    RHS = zeros(N*(nx+nu),(nx+nu)^2);
    for kk = 1:N
        tau_opt = [x_opt(:,kk);u_opt(:,kk)];
        tau_ref = [desired_state(:,ii+kk-1);
                   zeros(nu,1)];
        RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,:) = 0.5* (kron((-tau_opt+tau_ref)',eye(nx+nu)) + kron(eye(nx+nu),(-tau_opt+tau_ref)'));
    end
    % solve for the du/dQ11
    for kk = 1:N
        ocp.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,1),kk-1); % set the 1st column of RHS, which correspond to Q11
    end
    ocp.set('cost_y_ref_e', zeros(nx,1));
    
    ocp.solve;
    
    du_dQ11_new = ocp.get('u',0);
    
    % solve for the du/dQ22
    for kk = 1:N
        ocp.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,1*(nu+nx)+2),kk-1); % set the 6th column of RHS, which correspond to Q22
    end
    ocp.set('cost_y_ref_e', zeros(nx,1));
    
    ocp.solve;
    
    du_dQ22_new = ocp.get('u',0);
    
    % solve for the du/dQ33
    for kk = 1:N
        ocp.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,2*(nu+nx)+3),kk-1); % set the 15th column of RHS, which correspond to Q33
    end
    ocp.set('cost_y_ref_e', zeros(nx,1));
    
    ocp.solve;
    
    du_dQ33_new = ocp.get('u',0);
    
    
    % solve for the du/dR11
    for kk = 1:N
        ocp.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,3*(nu+nx)+4),kk-1); % set the 29th column of RHS, which correspond to R11
    end
    ocp.set('cost_y_ref_e', zeros(nx,1));
    
    ocp.solve;
    
    du_dR11_new = ocp.get('u',0);
    
    % solve for the du/dR22
    for kk = 1:N
        ocp.set('cost_y_ref',inv(W)*RHS((nx+nu)*(kk-1)+1:(nx+nu)*kk,4*(nu+nx)+5),kk-1); % set the 36th column of RHS, which correspond to R22
    end
    ocp.set('cost_y_ref_e', zeros(nx,1));
    
    ocp.solve;
    
    du_dR22_new = ocp.get('u',0);
    
    
    % only counting the diagonal elements
    du1_dQR = [du_dQ11_new(1) du_dQ22_new(1) du_dQ33_new(1) du_dR11_new(1) du_dR22_new(1)];
    du2_dQR = [du_dQ11_new(2) du_dQ22_new(2) du_dQ33_new(2) du_dR11_new(2) du_dR22_new(2)];
    
    du_dQR = [du1_dQR;du2_dQR];
end

function [du_dQR,timing_nu] = get_dQR_numerical(W,nx,nu,ocp_grad,x_sim,x_traj_init,u_traj_init,pi_traj_init,time_varying_y_ref,desired_state,ii,ocp_N)
    du_dQR=zeros(nu,nu+nx);
    time=[];
    for i = 1:nx
        % reset ocp parameter
        ocp_grad.set('constr_x0', x_sim(:,ii));
        ocp_grad.set('init_x', x_traj_init);
        ocp_grad.set('init_u', u_traj_init);
        ocp_grad.set('init_pi', pi_traj_init);
        
        for jj = 1:ocp_N
            ocp_grad.set('cost_y_ref', time_varying_y_ref(:,jj), jj-1);
        end
        time_varying_y_ref_N = desired_state(:,ii+ocp_N-1);
        ocp_grad.set('cost_y_ref_e', time_varying_y_ref_N);
        
        %perturb element of Q
        W_new1 = W;
        W_new1(i,i) = W_new1(i,i)*1.01;
        ocp_grad.set('cost_W',W_new1);
        tim1 = tic;
        ocp_grad.solve;
        time = [time toc(tim1)];
        u_new_p= vpa(ocp_grad.get('u',0));
        
        % reset ocp parameter
        ocp_grad.set('constr_x0', x_sim(:,ii));
        ocp_grad.set('init_x', x_traj_init);
        ocp_grad.set('init_u', u_traj_init);
        ocp_grad.set('init_pi', pi_traj_init);
        
        for jj = 1:ocp_N
            ocp_grad.set('cost_y_ref', time_varying_y_ref(:,jj), jj-1);
        end
        time_varying_y_ref_N = desired_state(:,ii+ocp_N-1);
        ocp_grad.set('cost_y_ref_e', time_varying_y_ref_N);
        W_new2 = W;
        W_new2(i,i) = W_new2(i,i)*0.99;
        ocp_grad.set('cost_W',W_new2);
        tim2 = tic;
        ocp_grad.solve;
        time = [time toc(tim2)];
        u_new_m= vpa(ocp_grad.get('u',0));
        W_new1(i,i)-W_new2(i,i);
        du_dQii_new= (u_new_p-u_new_m)/(W_new1(i,i)-W_new2(i,i));
        du_dQR(:,i) = du_dQii_new;
    end
    for i = 1:nu
        % reset ocp parameter
        ocp_grad.set('constr_x0', x_sim(:,ii));
        ocp_grad.set('init_x', x_traj_init);
        ocp_grad.set('init_u', u_traj_init);
        ocp_grad.set('init_pi', pi_traj_init);
        
        for jj = 1:ocp_N
            ocp_grad.set('cost_y_ref', time_varying_y_ref(:,jj), jj-1);
        end
        time_varying_y_ref_N = desired_state(:,ii+ocp_N-1);
        ocp_grad.set('cost_y_ref_e', time_varying_y_ref_N);

        %perturb element of R
        W_new1 = W;
        W_new1(nx+i,nx+i) = W_new1(nx+i,nx+i)*1.01;
        ocp_grad.set('cost_W',W_new1);
        tim3 = tic;
        ocp_grad.solve;
        time = [time toc(tim3)];
        u_new_p= vpa(ocp_grad.get('u',0));
        
        % reset ocp parameter
        ocp_grad.set('constr_x0', x_sim(:,ii));
        ocp_grad.set('init_x', x_traj_init);
        ocp_grad.set('init_u', u_traj_init);
        ocp_grad.set('init_pi', pi_traj_init);
        
        for jj = 1:ocp_N
            ocp_grad.set('cost_y_ref', time_varying_y_ref(:,jj), jj-1);
        end
        time_varying_y_ref_N = desired_state(:,ii+ocp_N-1);
        ocp_grad.set('cost_y_ref_e', time_varying_y_ref_N);

        W_new2 = W;
        W_new2(nx+i,nx+i) = W_new2(nx+i,nx+i)*0.99;
        ocp_grad.set('cost_W',W_new2);
        tim4 = tic;
        ocp_grad.solve;
        time = [time toc(tim4)];
        u_new_m= vpa(ocp_grad.get('u',0));
    
        du_dRii_new= (u_new_p-u_new_m)/(W_new1(nx+i,nx+i)-W_new2(nx+i,nx+i));
        du_dQR(:,nx+i) = du_dRii_new;
    end
    timing_nu = sum(time);
end
