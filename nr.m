% industrial_reactors_highres_full_fixed_complete.m
% Single-file production-grade reactor simulator (MATLAB/Octave)
% - Includes geometry setters, industrial params with defaults
% - Axisymmetric Laplacian (Kronecker), rhs_packed, analytic Jacobian
% - Circular HistoryCache, RK4 adaptive, SDIRK(3), Rosenbrock integrators
% - LU caching with ILU+GMRES fallback, memory-safe checkpoint wrapper
% Save as industrial_reactors_highres_full_fixed_complete.m and run.

##close all force; clc; clearvars -except; rand('state',1);
##
##function R = set_packed_bed_geometry(L, Rrad, Nz, Nr)
##R = struct();
##R.geometry.L = L; R.geometry.R = Rrad;
##R.grid.Nz = Nz; R.grid.Nr = Nr; R.Ncells = Nz * Nr;
##R.dx = L / (Nz - 1); R.dr = Rrad / (Nr - 1);
##R.ic.c0 = 0.5 * ones(R.Ncells,1); R.ic.T0 = 300.0 * ones(R.Ncells,1);
##R.delay.h = 0.3; R.history.grid = 0.005;
##R.bc.inlet.c = 0.8; R.bc.inlet.T = 300.0;
##end
##
##function R = set_multi_tubular_geometry(Ltube, Rtube, Nz, Nr)
##R = struct();
##R.geometry.tube_length = Ltube; R.geometry.tube_radius = Rtube;
##R.grid.Nz = Nz; R.grid.Nr = Nr; R.Ncells = Nz * Nr;
##R.dx = Ltube / (Nz - 1); R.dr = Rtube / (Nr - 1);
##R.ic.c0 = 0.4 * ones(R.Ncells,1); R.ic.T0 = 320.0 * ones(R.Ncells,1);
##R.delay.h = 0.15; R.history.grid = 0.005;
##R.bc.inlet.c = 0.6; R.bc.inlet.T = 320.0;
##end
##
##%% ---------------- Industrial params with defaults ----------------
##function s = default_solver_settings(method, atol, rtol, max_step, interp_grid, jacobian_reuse_tol)
##if nargin < 1, method = 'rosenbrock'; end
##if nargin < 2, atol = 1e-7; end
##if nargin < 3, rtol = 1e-7; end
##if nargin < 4, max_step = 0.02; end
##if nargin < 5, interp_grid = 0.005; end
##if nargin < 6, jacobian_reuse_tol = 1e-8; end
##s = struct(); s.method = method; s.opts = struct('atol',atol,'rtol',rtol,'max_step',max_step,'interp_grid',interp_grid); s.jacobian_reuse_tol = jacobian_reuse_tol;
##end
##
##function R = set_packed_bed_industrial_params(R, phys)
##if nargin < 2, phys = struct(); end
##phys.Dc = getfielddef(phys,'Dc',1.0e-5);
##phys.DT = getfielddef(phys,'DT',5.0e-6);
##phys.k1 = getfielddef(phys,'k1',10.0);
##phys.Ea = getfielddef(phys,'Ea',6.0e4);
##phys.k2 = getfielddef(phys,'k2',0.5);
##phys.beta = getfielddef(phys,'beta',5.0);
##phys.gamma = getfielddef(phys,'gamma',1.0);
##phys.porosity = getfielddef(phys,'porosity',0.4);
##phys.rho = getfielddef(phys,'rho',1.2);
##phys.cp = getfielddef(phys,'cp',1000.0);
##R.phys = phys;
##end
##
##function R = set_multi_tubular_industrial_params(R, phys)
##if nargin < 2, phys = struct(); end
##phys.Dc = getfielddef(phys,'Dc',5.0e-6);
##phys.DT = getfielddef(phys,'DT',2.0e-6);
##phys.k1 = getfielddef(phys,'k1',1.2);
##phys.Ea = getfielddef(phys,'Ea',4.5e4);
##phys.k2 = getfielddef(phys,'k2',0.15);
##phys.beta = getfielddef(phys,'beta',2.0);
##phys.gamma = getfielddef(phys,'gamma',0.6);
##phys.porosity = getfielddef(phys,'porosity',0.45);
##phys.rho = getfielddef(phys,'rho',1.0);
##phys.cp = getfielddef(phys,'cp',1200.0);
##R.phys = phys;
##end
##
##%% ---------------- Axisymmetric Laplacian (Kronecker) ----------------
##function L = build_axisymmetric_laplacian_kron(Nz, Nr, dz, dr)
##% faster Kronecker construction
##ex = spdiags([ones(Nz,1), -2*ones(Nz,1), ones(Nz,1)], [-1,0,1], Nz, Nz);
##ey = spdiags([ones(Nr,1), -2*ones(Nr,1), ones(Nr,1)], [-1,0,1], Nr, Nr);
##Ix = speye(Nz); Iy = speye(Nr);
##Lz = ex / (dz*dz); Lr = ey / (dr*dr);
##L = kron(Iy, Lz) + kron(Lr, Ix);
##end
##
##%% ---------------- Reactor RHS (packed) ----------------
##function dX = rhs_packed(t, X, X_tau, Lmat, P)
##% X = [c; T], X_tau same shape
##N = P.Ncells;
##c = X(1:N); T = X(N+1:end);
##% Arrhenius-like rate (lumped)
##Rrate = P.phys.k1 .* c .* exp(-(P.phys.Ea/8.314) .* (1./T - 1/300));
##dc = P.phys.Dc * (Lmat * c) - Rrate - P.phys.k2 .* c;
##dT = P.phys.DT * (Lmat * T) + P.phys.beta .* Rrate - P.phys.gamma .* (X_tau(N+1:end) - P.bc.inlet.T);
##dX = [dc; dT];
##end
##
##%% ---------------- Analytic Jacobian (sparse block) ----------------
##function J = analytic_jacobian_reactor(x, params, Lmat)
##N = params.Ncells;
##c = x(1:N); T = x(N+1:end);
##k1 = params.phys.k1; k2 = params.phys.k2; beta = params.phys.beta; gamma = params.phys.gamma;
##expT = exp(-(params.phys.Ea/8.314) .* (1./T - 1/300));
##dr_dc = k1 .* expT;
##dexp_dT = expT .* (params.phys.Ea/8.314) .* (1./(T.^2));
##dr_dT = k1 .* c .* dexp_dT;
##J11 = params.phys.Dc * Lmat - spdiags(k2 * ones(N,1) + dr_dc, 0, N, N);
##J12 = - spdiags(dr_dT, 0, N, N);
##J21 = spdiags(beta * dr_dc, 0, N, N);
##J22 = params.phys.DT * Lmat - spdiags(gamma * ones(N,1), 0, N, N);
##J = [J11, J12; J21, J22];
##end
##
##%% ---------------- HistoryCache (circular) ----------------
##function H = HistoryCache_create(t_init, x_init, maxlen)
##if nargin < 3, maxlen = 200000; end
##nt = numel(t_init); nstate = size(x_init,2);
##H.maxlen = max(maxlen, nt);
##H.tbuf = nan(H.maxlen,1); H.xbuf = nan(H.maxlen, nstate);
##H.head = 1; H.tail = nt; H.count = nt;
##H.tbuf(1:nt) = t_init(:); H.xbuf(1:nt,:) = x_init;
##end
##
##function H = HistoryCache_append(H, t_new, x_new)
##if H.count < H.maxlen
##  H.tail = H.tail + 1; H.tbuf(H.tail) = t_new; H.xbuf(H.tail,:) = x_new(:).'; H.count = H.count + 1;
##else
##  H.head = H.head + 1; if H.head > H.maxlen, H.head = 1; end
##  H.tail = H.head + H.count - 1; if H.tail > H.maxlen, H.tail = H.tail - H.maxlen; end
##  H.tbuf(H.tail) = t_new; H.xbuf(H.tail,:) = x_new(:).';
##end
##end
##
##function xq = HistoryCache_eval(H, tq)
##if H.count == 0, error('HistoryCache is empty'); end
##if H.head <= H.tail
##  tvec = H.tbuf(H.head:H.tail); xmat = H.xbuf(H.head:H.tail,:);
##else
##  tvec = [H.tbuf(H.head:end); H.tbuf(1:H.tail)]; xmat = [H.xbuf(H.head:end,:); H.xbuf(1:H.tail,:)];
##end
##xq = interp1(tvec, xmat, tq, 'pchip')';
##end
##
##%% ---------------- run_and_checkpoint wrapper ----------------
##function [Tchk,Xchk,final_state,stats] = run_and_checkpoint(solver_fn, f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, label, save_every, max_nstate_for_lu, history_margin)
##
####
####
####if nargin < 11, history_margin = 10; end
####if nargin < 10, max_nstate_for_lu = 5e4; end
####if nargin < 9, save_every = 50; end
####
####est_steps = ceil((tspan(2)-tspan(1)) / max(1e-6, opts.max_step)) + 10;
####est_store = max(10, ceil(est_steps / max(1, save_every)));
####Tchk = zeros(est_store,1); Xchk = zeros(est_store, n_state); idx_store = 0;
####
####interp_grid = getfielddef(opts,'interp_grid',0.01);
####hist_max = getfielddef(opts,'history_maxlen', ceil(h_delay / interp_grid) + history_margin);
####opts.history_maxlen = hist_max;
####
####if n_state > max_nstate_for_lu, opts.force_ilu = true; else opts.force_ilu = getfielddef(opts,'force_ilu', false); end
####
####step_counter = 0;
####
####
####[Tfull,Xfull,stats] = solver_fn(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, @accept_cb);
####
####final_state = Xfull(end,:).';
####if idx_store == 0
####  Tchk = Tfull(end); Xchk = Xfull(end,:);
####else
####  Tchk = Tchk(1:idx_store); Xchk = Xchk(1:idx_store,:);
####end
####
####clear Tfull Xfull; try pack; catch, end
####fprintf('%s done: stored %d checkpoints\n', label, numel(Tchk));
####end
####
####function accept_cb(t_step, x_step,Tc)
####  step_counter = step_counter + 1;
####  if mod(step_counter, save_every) == 0 || t_step >= tspan(2)-1e-12
####    idx_store = idx_store + 1;
####    if idx_store > size(Tchk,1)
####      grow = max(ceil(est_store/2),10);
####      Tchk = [Tchk; zeros(grow,1)]; Xchk = [Xchk; zeros(grow, n_state)];
####    end
####    Tchk(idx_store) = t_step; Xchk(idx_store,:) = x_step(:).';
####  end
##
## % RUN_AND_CHECKPOINT  Run a solver with downsampling via external callback
##    %
##    % Inputs:
##    %   solver_fn   - function handle to integrator (e.g. @sdirk)
##    %   f_rhs       - RHS of DDE
##    %   x0_fun      - history function
##    %   h_delay     - delay
##    %   tspan       - [t0 tfinal]
##    %   n_state     - number of state variables
##    %   extra_args  - cell array of extra arguments to f_rhs
##    %   opts        - solver options structure
##    %   label       - string for progress message
##    %   save_every  - store every 'save_every' steps (default 50)
##    %
##    % Outputs:
##    %   Tchk        - times of stored checkpoints
##    %   Xchk        - states at those times
##    %   final_state - final state vector
##    %   stats       - statistics structure from solver
##
##    if nargin < 10
##        save_every = 50;
##    end
##
##    % Estimate maximum number of steps and preallocate
##    est_max_steps = ceil((tspan(2)-tspan(1)) / max(1e-6, opts.max_step)) + 10;
##    est_store = ceil(est_max_steps / save_every) + 2;
##
##    % Global state for callback
##    global callback_state
##    callback_state = struct();
##    callback_state.Tchk = zeros(est_store, 1);
##    callback_state.Xchk = zeros(est_store, n_state);
##    callback_state.idx_store = 0;
##    callback_state.step_counter = 0;
##    callback_state.save_every = save_every;
##    callback_state.tspan = tspan;
##    callback_state.est_store = est_store;
##    callback_state.n_state = n_state;
##
##    % Call solver with external callback
##    [Tfull, Xfull, stats] = solver_fn(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, @run_and_checkpoint_callback);
##
##    % Extract final state
##    final_state = Xfull(end, :).';
##
##    % Retrieve checkpoint arrays from global state
##    idx_store = callback_state.idx_store;
##    if idx_store == 0
##        Tchk = Tfull(end);
##        Xchk = Xfull(end, :);
##    else
##        Tchk = callback_state.Tchk(1:idx_store);
##        Xchk = callback_state.Xchk(1:idx_store, :);
##    end
##
##    % Clear global to avoid interference with next run
##    clear global callback_state
##
##    % Free memory
##    clear Tfull Xfull;
##    try pack; catch; end
##
##    fprintf('%s done: stored %d checkpoints\n', label, numel(Tchk));
##end
##
##% ----------------------------------------------------------------------
##% External callback function used by run_and_checkpoint
##% ----------------------------------------------------------------------
##function run_and_checkpoint_callback(t_step, x_step)
##    global callback_state
##
##    callback_state.step_counter = callback_state.step_counter + 1;
##    if mod(callback_state.step_counter, callback_state.save_every) == 0 || ...
##       t_step >= callback_state.tspan(2) - 1e-12
##        callback_state.idx_store = callback_state.idx_store + 1;
##        idx = callback_state.idx_store;
##
##        % Grow arrays if necessary
##        if idx > size(callback_state.Tchk, 1)
##            extra = callback_state.est_store;
##            callback_state.Tchk = [callback_state.Tchk; zeros(extra, 1)];
##            callback_state.Xchk = [callback_state.Xchk; zeros(extra, callback_state.n_state)];
##        end
##
##        callback_state.Tchk(idx) = t_step;
##        callback_state.Xchk(idx, :) = x_step(:).';
##    end
##end
##
##
##
##%% ---------------- RK4 adaptive (callback supported) ----------------
##function [T_out,X_out,stats] = simulate_dde_rk4_adaptive(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
##if nargin < 8, accept_callback = []; end
##at = getfielddef(opts,'atol',1e-6); rt = getfielddef(opts,'rtol',1e-5); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
##interp_grid = getfielddef(opts,'interp_grid',0.05);
##h_eff = max(h_delay,1e-8); hist_times = (-h_eff:interp_grid:0); hist_states = zeros(numel(hist_times), n_state);
##for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
##H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
##t = tspan(1); x = x0_fun(0); x = x(:);
##T_out = t; X_out = x.'; stats.rhs_evals = 0; stats.steps = 0;
##h_try = min(max_step, 0.01*(tspan(2)-tspan(1)) + 1e-2);
##while t < tspan(2)-1e-12
##  h_try = min(h_try, tspan(2)-t);
##  x_tau = HistoryCache_eval(H, t - h_delay);
##  k1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
##  x2 = x + 0.5*h_try*k1; x_tau2 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
##  k2 = f_rhs(t + 0.5*h_try, x2, x_tau2); stats.rhs_evals = stats.rhs_evals + 1;
##  x3 = x + 0.5*h_try*k2; x_tau3 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
##  k3 = f_rhs(t + 0.5*h_try, x3, x_tau3); stats.rhs_evals = stats.rhs_evals + 1;
##  x4 = x + h_try*k3; x_tau4 = HistoryCache_eval(H, t + h_try - h_delay);
##  k4 = f_rhs(t + h_try, x4, x_tau4); stats.rhs_evals = stats.rhs_evals + 1;
##  x_rk4 = x + (h_try/6)*(k1 + 2*k2 + 2*k3 + k4);
##  x_heun = x + (h_try/2)*(k1 + k4);
##  err_est = norm(x_rk4 - x_heun, inf);
##  tol = at + rt * max(norm(x,inf), norm(x_rk4,inf));
##  if err_est <= tol
##    t = t + h_try; x = x_rk4;
##    T_out(end+1,1) = t; X_out(end+1,:) = x.'; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
##    if ~isempty(accept_callback), accept_callback(t,x); end
##  end
##  if err_est == 0, fac = 2.0; else fac = 0.8 * (tol/err_est)^0.2; end
##  h_try = min(max_step, max(1e-8, fac * h_try));
##end
##end
##
##%% ---------------- SDIRK(3) (callback supported) ----------------
##function [T_out,X_out,stats] = simulate_dde_sdirk(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
##if nargin < 8, accept_callback = []; end
##at = getfielddef(opts,'atol',1e-7); rt = getfielddef(opts,'rtol',1e-7); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
##interp_grid = getfielddef(opts,'interp_grid',0.02); Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
##gamma = 0.435866521508459; A = [gamma 0 0; 0.5-gamma gamma 0; 2*gamma 1-4*gamma gamma]; b = [2*gamma;1-4*gamma;2*gamma];
##h_eff = max(h_delay,1e-8); hist_times = (-h_eff:interp_grid:0); hist_states = zeros(numel(hist_times), n_state);
##for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
##H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
##t = tspan(1); x = x0_fun(0); x = x(:); T_out = t; X_out = x.'; stats.rhs_evals = 0; stats.newton_iters = 0; stats.linear_solves = 0; stats.steps = 0;
##Lmat = extra_args{1}; params = extra_args{2};
##J_prev = []; LU = struct('L',[],'U',[],'P',[],'Q',[],'R',[],'ilu',false);
##h = min(max_step, 1e-2); force_ilu = getfielddef(opts,'force_ilu',false);
##while t < tspan(2)-1e-12
##  h = min(h, tspan(2)-t); K = zeros(n_state,3);
##  for stage = 1:3
##    xi = x;
##    for newt = 1:8
##      ti = t + sum(A(stage,1:stage))*h; x_tau = HistoryCache_eval(H, ti - h_delay);
##      fval = f_rhs(ti, xi, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
##      R = xi - x - h*(A(stage,stage)*fval + sum(A(stage,1:stage-1).*K(:,1:stage-1),2));
##      J = analytic_jacobian_reactor(xi, params, Lmat); M = speye(n_state) - h*A(stage,stage)*J;
##      recompute = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
##      if recompute
##        if force_ilu
##          try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.P=[]; LU.Q=[]; LU.R=[]; LU.ilu=true;
##          catch, [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false; end
##        else
##          try [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false;
##          catch, try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.P=[]; LU.Q=[]; LU.R=[]; LU.ilu=true; catch LU.ilu=false; end
##          end
##        end
##        J_prev = J;
##      end
##      RHS = -R;
##      if LU.ilu
##        try [delta,flag] = gmres(M, RHS, [], 1e-8, 200, LU.L, LU.U); if flag ~= 0, delta = M \ RHS; end
##        catch delta = M \ RHS; end
##      else
##        if ~isempty(LU.R), y = LU.P * RHS; z = LU.L \ y; w = LU.U \ z; w = LU.R \ w; delta = LU.Q * w;
##        else y = LU.P * RHS; z = LU.L \ y; w = LU.U \ z; delta = LU.Q * w; end
##      end
##      xi = xi + delta; stats.newton_iters = stats.newton_iters + 1; stats.linear_solves = stats.linear_solves + 1;
##      if norm(delta, inf) < 1e-8, break; end
##    end
##    K(:,stage) = fval;
##  end
##  x_new = x + h*(K * b);
##  err_est = norm(x_new - X_out(end,:).', inf);
##  tol = at + rt * max(norm(x,inf), norm(x_new,inf));
##  if err_est <= tol
##    t = t + h; x = x_new; T_out(end+1,1) = t; X_out(end+1,:) = x.'; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
##    if ~isempty(accept_callback), accept_callback(t,x); end
##    if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
##    h = min(max_step, max(1e-8, fac*h));
##  else
##    h = max(1e-8, 0.5*h);
##  end
##end
##end
##
##%% ---------------- Rosenbrock (callback supported) ----------------
##function [T_out,X_out,stats] = simulate_dde_rosenbrock(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
##if nargin < 8, accept_callback = []; end
##at = getfielddef(opts,'atol',1e-7); rt = getfielddef(opts,'rtol',1e-7); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
##interp_grid = getfielddef(opts,'interp_grid',0.02); Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
##gamma = 0.5;
##h_eff = max(h_delay,1e-8); hist_times = (-h_eff:interp_grid:0); hist_states = zeros(numel(hist_times), n_state);
##for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
##H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
##t = tspan(1); x = x0_fun(0); x = x(:); T_out = t; X_out = x.'; stats.rhs_evals = 0; stats.linear_solves = 0; stats.steps = 0;
##Lmat = extra_args{1}; params = extra_args{2};
##J_prev = []; LU = struct('L',[],'U',[],'P',[],'Q',[],'R',[],'ilu',false);
##h = min(max_step, 1e-2); force_ilu = getfielddef(opts,'force_ilu',false);
##while t < tspan(2)-1e-12
##  h = min(h, tspan(2)-t);
##  x_tau = HistoryCache_eval(H, t - h_delay);
##  J = analytic_jacobian_reactor(x, params, Lmat); M = speye(n_state) - gamma*h*J;
##  recompute = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
##  if recompute
##    if force_ilu
##      try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.ilu=true; catch, [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false; end
##    else
##      try [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false;
##      catch, try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.ilu=true; catch LU.ilu=false; end
##      end
##    end
##    J_prev = J;
##  end
##  f1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
##  if LU.ilu, [k1,fl]=gmres(M,f1,[],1e-8,200,LU.L,LU.U); if fl~=0, k1=M\f1; end
##  else y=LU.P*f1; z=LU.L\y; w=LU.U\z; if ~isempty(LU.R), w=LU.R\w; end; k1=LU.Q*w; end
##  t2 = t + 0.5*h; x2 = x + 0.5*h*k1; x2_tau = HistoryCache_eval(H, t2 - h_delay);
##  f2 = f_rhs(t2, x2, x2_tau); stats.rhs_evals = stats.rhs_evals + 1;
##  rhs2 = f2 + (1/(2*h))*k1;
##  if LU.ilu, [k2,fl]=gmres(M,rhs2,[],1e-8,200,LU.L,LU.U); if fl~=0, k2=M\rhs2; end
##  else y2=LU.P*rhs2; z2=LU.L\y2; w2=LU.U\z2; if ~isempty(LU.R), w2=LU.R\w2; end; k2=LU.Q*w2; end
##  x_new = x + h*(k1 + k2)/2;
##  err_est = norm(x_new - (x + h*k1), inf);
##  tol = at + rt * max(norm(x,inf), norm(x_new,inf));
##  if err_est <= tol
##    t = t + h; x = x_new; T_out(end+1,1) = t; X_out(end+1,:) = x.'; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
##    if ~isempty(accept_callback), accept_callback(t,x); end
##    if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
##    h = min(max_step, max(1e-8, fac*h));
##  else
##    h = max(1e-8, 0.5*h);
##  end
##end
##end
##
##%% ---------------- Plotting from checkpoints ----------------
##function create_3d_plots_from_checkpoints(P, Lmat, Tchk, Xchk, title_prefix)
##if isempty(Tchk) || isempty(Xchk), warning('No checkpoints to plot for %s', title_prefix); return; end
##Nz = P.grid.Nz; Nr = P.grid.Nr; N = P.Ncells;
##cf = Xchk(end,1:N)'; tf = Xchk(end,N+1:end)';
##C = reshape(cf,[Nr Nz]); Tm = reshape(tf,[Nr Nz]);
##z = linspace(0, getfielddef(P,'geometry.L',getfielddef(P,'geometry.tube_length',1.0)), Nz);
##r = linspace(0, getfielddef(P,'geometry.R',getfielddef(P,'geometry.tube_radius',1.0)), Nr);
##max_cells_for_isosurface = 200000;
##if P.Ncells <= max_cells_for_isosurface
##  nang = 36; theta = linspace(0,2*pi,nang);
##  Xcart = zeros(Nr,nang,Nz); Ycart = Xcart; Zcart = zeros(Nr,nang,Nz);
##  Ccart = zeros(Nr,nang,Nz); Tcart = zeros(Nr,nang,Nz);
##  for iz=1:Nz, for ir=1:Nr, for it=1:nang
##    Xcart(ir,it,iz)=r(ir)*cos(theta(it)); Ycart(ir,it,iz)=r(ir)*sin(theta(it)); Zcart(ir,it,iz)=z(iz);
##    Ccart(ir,it,iz)=C(ir,iz); Tcart(ir,it,iz)=Tm(ir,iz);
##  end,end,end
##  figure('Name',[title_prefix ' Temp iso']); p=patch(isosurface(Xcart,Ycart,Zcart,Tcart,mean(Tcart(:)))); isonormals(Xcart,Ycart,Zcart,Tcart,p);
##  set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7); camlight; lighting gouraud; axis equal;
##else
##  figure('Name',[title_prefix ' Conc slice']); imagesc(z,r,C); colorbar; axis xy; xlabel('z (m)'); ylabel('r (m)');
##  figure('Name',[title_prefix ' Centerline T']); plot(z, Tm(ceil(Nr/2),:)); xlabel('z (m)'); ylabel('T (K)');
##end
##end
##
##%% ---------------- Utility ----------------
##function v = getfielddef(s,f,d)
##if isempty(s) || ~isstruct(s) || ~isfield(s,f), v = d; else v = s.(f); end
##end
##
##
##function industrial_reactors_highres_full_fixed_complete()
##
## try pack; catch, end
##
##% ---------------- User geometry and industrial parameters ----------------
##A = set_packed_bed_geometry(12.0,1.6,240,60);
##B = set_multi_tubular_geometry(18.0,0.15,400,32);
##
##A = set_packed_bed_industrial_params(A, struct('Dc',1e-5,'DT',5e-6,'k1',10,'Ea',6e4,'k2',0.5,'beta',5,'gamma',1,'porosity',0.4));
##B = set_multi_tubular_industrial_params(B, struct('Dc',5e-6,'DT',2e-6,'k1',1.2,'Ea',4.5e4,'k2',0.15,'beta',2,'gamma',0.6,'porosity',0.45));
##
##A.Tfinal = 60; B.Tfinal = 120;
##A.solver = default_solver_settings('rosenbrock',1e-7,1e-7,0.02,0.005,1e-8);
##B.solver = default_solver_settings('sdirk',1e-7,1e-7,0.02,0.005,1e-8);
##
##% ---------------- Spatial operators (reuse) ----------------
##L_A = build_axisymmetric_laplacian_kron(A.grid.Nz, A.grid.Nr, A.dx, A.dr);
##L_B = build_axisymmetric_laplacian_kron(B.grid.Nz, B.grid.Nr, B.dx, B.dr);
##
##x0A_fun = @(t) [A.ic.c0; A.ic.T0];
##x0B_fun = @(t) [B.ic.c0; B.ic.T0];
##
##% ---------------- Memory-safe run parameters ----------------
##save_every = 50;            % store every N accepted steps
##max_nstate_for_lu = 5e4;    % threshold to prefer ILU+GMRES
##history_margin = 10;        % extra history slots
##
##fprintf('Running high-resolution simulations (memory-safe mode)...\n');
##
##
###Fuel types:
####
####CNG
####LNG
####Propane
####Hydrogene
##
##
##% Reactor A: SDIRK reference (coarse tolerances to reduce cost)
##opts_refA = struct('atol',1e-8,'rtol',1e-8,'max_step',A.Tfinal/500,'interp_grid',A.history.grid,'jacobian_reuse_tol',A.solver.jacobian_reuse_tol);
##[TA_ref_chk, XA_ref_chk, finalA_ref, stats_refA] = run_and_checkpoint(@simulate_dde_sdirk, @(t,x,xt) rhs_packed(t,x,xt,L_A,A), x0A_fun, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {L_A,A}, opts_refA, 'SDIRK-A', save_every, max_nstate_for_lu, history_margin);
##
##% Reactor A: Rosenbrock
##opts_rosA = A.solver.opts; opts_rosA.jacobian_reuse_tol = A.solver.jacobian_reuse_tol;
##[TA_ros_chk, XA_ros_chk, finalA_ros, stats_rosA] = run_and_checkpoint(@simulate_dde_rosenbrock, @(t,x,xt) rhs_packed(t,x,xt,L_A,A), x0A_fun, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {L_A,A}, opts_rosA, 'ROSENBROCK-A', save_every, max_nstate_for_lu, history_margin);
##
##clear TA_ref_chk XA_ref_chk; try pack; catch, end
##
##% Reactor B: SDIRK
##opts_sdirkB = B.solver.opts; opts_sdirkB.jacobian_reuse_tol = B.solver.jacobian_reuse_tol;
##[TB_sdirk_chk, XB_sdirk_chk, finalB_sdirk, stats_sdirkB] = run_and_checkpoint(@simulate_dde_sdirk, @(t,x,xt) rhs_packed(t,x,xt,L_B,B), x0B_fun, B.delay.h, [0 B.Tfinal], 2*B.Ncells, {L_B,B}, opts_sdirkB, 'SDIRK-B', save_every, max_nstate_for_lu, history_margin);
##
##% Save compact results
##save('results_A_rosenbrock.mat','TA_ros_chk','XA_ros_chk','finalA_ros','-v7.3');
##save('results_B_sdirk.mat','TB_sdirk_chk','XB_sdirk_chk','finalB_sdirk','-v7.3');
##
##% Plot from checkpoints
##try
##  create_3d_plots_from_checkpoints(A, L_A, TA_ros_chk, XA_ros_chk, 'Packed-bed reactor (A)');
##  create_3d_plots_from_checkpoints(B, L_B, TB_sdirk_chk, XB_sdirk_chk, 'Multi-tubular reactor (B)');
##catch ME
##  warning('Plotting failed: %s', ME.message);
##end
##
##fprintf('All runs complete; checkpoints saved.\n');
##end
##
##%% ---------------- Geometry setters ----------------
##
##industrial_reactors_highres_full_fixed_complete()
##% End of file

##  https://contra.com/mr_pavl_mazniker_edxz30gh?referralExperimentNid=DEFAULT_REFERRAL_PROGRAM&referrerUsername=mr_pavl_mazniker_edxz30gh
##  +380990535261
##  https://diag.net/u/u6r3ondjie0w0l8138bafm095b
##  https://github.com/goodengineer
##  https://gust.com/companies/mr-pavel-s-startup
##  https://orcid.org/0000-0001-8184-8166
##    https://willwork781147312.wordpress.com/portfolio/cp/
##    https://www.mathworks.com/matlabcentral/profile/authors/1643504
##    https://www.researchgate.net/profile/Pavel-Mazniker
##    https://nanohub.org/members/130066
##    https://pavelmaz.blogspot.com/
##    ko-fi.com/mrpaul61277

##    https://independant.academia.edu/PavelMazniker
##    https://www.youtube.com/channel/UCC__7jMOAHak0MVkUFtmO-w
##    https://scholar.google.co.uk/citations?user=cnTX270AAAAJ&hl=en
##    https://flippa.com/users/4752309



% industrial_reactors_highres_memory_safe.m
% Single-file, memory-safe reactor simulator (MATLAB / Octave)
% - Main entry: industrial_reactors_highres_memory_safe()
% - Memory-safety: integrators stream accepted steps to a callback; they do NOT retain full trajectories
% - Checkpoint wrapper collects downsampled checkpoints and final state only
% - Sparse Kronecker Laplacian, analytic Jacobian, LU caching with ILU+GMRES fallback
% - Use this file as a drop-in: save and run industrial_reactors_highres_memory_safe()


##%% ------------------------------------------------------------------------
##%% run_and_checkpoint: wrapper that collects downsampled checkpoints only
##%% ------------------------------------------------------------------------
##function [Tchk,Xchk,final_state,stats] = run_and_checkpoint(solver_fn, f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, label, save_every, max_nstate_for_lu, history_margin)
##
##if nargin < 11, history_margin = 10; end
##if nargin < 10, max_nstate_for_lu = 5e4; end
##if nargin < 9, save_every = 50; end
##
##% conservative estimate for preallocation
##est_steps = ceil((tspan(2)-tspan(1)) / max(1e-6, getfielddef(opts,'max_step', (tspan(2)-tspan(1))/100))) + 10;
##est_store = max(10, ceil(est_steps / max(1, save_every)));
##
##% prepare callback state (local struct, not global)
##cb.Tchk = zeros(est_store,1);
##cb.Xchk = zeros(est_store, n_state);
##cb.idx_store = 0;
##cb.step_counter = 0;
##cb.save_every = save_every;
##cb.tspan = tspan;
##cb.est_store = est_store;
##cb.n_state = n_state;
##
##% adjust history buffer size
##interp_grid = getfielddef(opts,'interp_grid',0.01);
##opts.history_maxlen = getfielddef(opts,'history_maxlen', ceil(h_delay / interp_grid) + history_margin);
##
##% force ILU for very large systems
##if n_state > max_nstate_for_lu
##  opts.force_ilu = true;
##else
##  opts.force_ilu = getfielddef(opts,'force_ilu', false);
##end
##
##% define callback closure
##function accept_cb(t_step, x_step)
##  cb.step_counter = cb.step_counter + 1;
##  if mod(cb.step_counter, cb.save_every) == 0 || t_step >= cb.tspan(2) - 1e-12
##    cb.idx_store = cb.idx_store + 1;
##    if cb.idx_store > size(cb.Tchk,1)
##      extra = cb.est_store
##      cb.Tchk = [cb.Tchk; zeros(extra,1)];
##      cb.Xchk = [cb.Xchk; zeros(extra, cb.n_state)];
##    end
##    cb.Tchk(cb.idx_store) = t_step;
##    cb.Xchk(cb.idx_store, :) = x_step(:).';
##  end
##end
##
##% call solver: integrators return final_state and stats (they stream via callback)
##[final_state, stats] = solver_fn(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, @accept_cb);
##
##% extract checkpoints
##if cb.idx_store == 0
##  Tchk = final_state(1); % fallback: single time
##  Xchk = final_state(:).';
##else
##  Tchk = cb.Tchk(1:cb.idx_store);
##  Xchk = cb.Xchk(1:cb.idx_store, :);
##end
##
##fprintf('%s done: stored %d checkpoints\n', label, numel(Tchk));
##end
##
##%% ------------------------------------------------------------------------
##%% Geometry setters
##%% ------------------------------------------------------------------------
##function R = set_packed_bed_geometry(L, Rrad, Nz, Nr)
##R = struct();
##R.geometry.L = L; R.geometry.R = Rrad;
##R.grid.Nz = Nz; R.grid.Nr = Nr; R.Ncells = Nz * Nr;
##R.dx = L / (Nz - 1); R.dr = Rrad / (Nr - 1);
##R.ic.c0 = 0.5 * ones(R.Ncells,1); R.ic.T0 = 300.0 * ones(R.Ncells,1);
##R.delay.h = 0.3; R.history.grid = 0.005;
##R.bc.inlet.c = 0.8; R.bc.inlet.T = 300.0;
##end
##
##function R = set_multi_tubular_geometry(Ltube, Rtube, Nz, Nr)
##R = struct();
##R.geometry.tube_length = Ltube; R.geometry.tube_radius = Rtube;
##R.grid.Nz = Nz; R.grid.Nr = Nr; R.Ncells = Nz * Nr;
##R.dx = Ltube / (Nz - 1); R.dr = Rtube / (Nr - 1);
##R.ic.c0 = 0.4 * ones(R.Ncells,1); R.ic.T0 = 320.0 * ones(R.Ncells,1);
##R.delay.h = 0.15; R.history.grid = 0.005;
##R.bc.inlet.c = 0.6; R.bc.inlet.T = 320.0;
##end
##
##%% ------------------------------------------------------------------------
##%% Industrial params with defaults
##%% ------------------------------------------------------------------------
##function s = default_solver_settings(method, atol, rtol, max_step, interp_grid, jacobian_reuse_tol)
##if nargin < 1, method = 'rosenbrock'; end
##if nargin < 2, atol = 1e-7; end
##if nargin < 3, rtol = 1e-7; end
##if nargin < 4, max_step = 0.02; end
##if nargin < 5, interp_grid = 0.005; end
##if nargin < 6, jacobian_reuse_tol = 1e-8; end
##s = struct(); s.method = method; s.opts = struct('atol',atol,'rtol',rtol,'max_step',max_step,'interp_grid',interp_grid); s.jacobian_reuse_tol = jacobian_reuse_tol;
##end
##
##function R = set_packed_bed_industrial_params(R, phys)
##if nargin < 2, phys = struct(); end
##phys.Dc = getfielddef(phys,'Dc',1.0e-5);
##phys.DT = getfielddef(phys,'DT',5.0e-6);
##phys.k1 = getfielddef(phys,'k1',10.0);
##phys.Ea = getfielddef(phys,'Ea',6.0e4);
##phys.k2 = getfielddef(phys,'k2',0.5);
##phys.beta = getfielddef(phys,'beta',5.0);
##phys.gamma = getfielddef(phys,'gamma',1.0);
##phys.porosity = getfielddef(phys,'porosity',0.4);
##phys.rho = getfielddef(phys,'rho',1.2);
##phys.cp = getfielddef(phys,'cp',1000.0);
##R.phys = phys;
##end
##
##function R = set_multi_tubular_industrial_params(R, phys)
##if nargin < 2, phys = struct(); end
##phys.Dc = getfielddef(phys,'Dc',5.0e-6);
##phys.DT = getfielddef(phys,'DT',2.0e-6);
##phys.k1 = getfielddef(phys,'k1',1.2);
##phys.Ea = getfielddef(phys,'Ea',4.5e4);
##phys.k2 = getfielddef(phys,'k2',0.15);
##phys.beta = getfielddef(phys,'beta',2.0);
##phys.gamma = getfielddef(phys,'gamma',0.6);
##phys.porosity = getfielddef(phys,'porosity',0.45);
##phys.rho = getfielddef(phys,'rho',1.0);
##phys.cp = getfielddef(phys,'cp',1200.0);
##R.phys = phys;
##end
##
##%% ------------------------------------------------------------------------
##%% Axisymmetric Laplacian (Kronecker)
##%% ------------------------------------------------------------------------
##function L = build_axisymmetric_laplacian_kron(Nz, Nr, dz, dr)
##ex = spdiags([ones(Nz,1), -2*ones(Nz,1), ones(Nz,1)], [-1,0,1], Nz, Nz);
##ey = spdiags([ones(Nr,1), -2*ones(Nr,1), ones(Nr,1)], [-1,0,1], Nr, Nr);
##Ix = speye(Nz); Iy = speye(Nr);
##Lz = ex / (dz*dz); Lr = ey / (dr*dr);
##L = kron(Iy, Lz) + kron(Lr, Ix)
##end
##
##%% ------------------------------------------------------------------------
##%% Reactor RHS (packed)
##%% ------------------------------------------------------------------------
##function dX = rhs_packed(t, X, X_tau, Lmat, P)
##N = P.Ncells;
##c = X(1:N); T = X(N+1:end);
##Rrate = P.phys.k1 .* c .* exp(-(P.phys.Ea/8.314) .* (1./T - 1/300));
##dc = P.phys.Dc * (Lmat * c) - Rrate - P.phys.k2 .* c;
##dT = P.phys.DT * (Lmat * T) + P.phys.beta .* Rrate - P.phys.gamma .* (X_tau(N+1:end) - P.bc.inlet.T);
##dX = [dc; dT];
##end
##
##%% ------------------------------------------------------------------------
##%% Analytic Jacobian (sparse block)
##%% ------------------------------------------------------------------------
##function J = analytic_jacobian_reactor(x, params, Lmat)
##N = params.Ncells
##c = x(1:N); T = x(N+1:end);
##k1 = params.phys.k1; k2 = params.phys.k2; beta = params.phys.beta; gamma = params.phys.gamma;
##expT = exp(-(params.phys.Ea/8.314) .* (1./T - 1/300));
##dr_dc = k1 .* expT;
##dexp_dT = expT .* (params.phys.Ea/8.314) .* (1./(T.^2));
##dr_dT = k1 .* c .* dexp_dT;
##J11 = params.phys.Dc * Lmat - spdiags(k2 * ones(N,1) + dr_dc, 0, N, N);
##J12 = - spdiags(dr_dT, 0, N, N);
##J21 = spdiags(beta * dr_dc, 0, N, N);
##J22 = params.phys.DT * Lmat - spdiags(gamma * ones(N,1), 0, N, N);
##J = [J11, J12; J21, J22];
##end
##
##%% ------------------------------------------------------------------------
##%% HistoryCache (circular)
##%% ------------------------------------------------------------------------
##function H = HistoryCache_create(t_init, x_init, maxlen)
##if nargin < 3, maxlen = 200000; end
##nt = numel(t_init); nstate = size(x_init,2);
##H.maxlen = max(maxlen, nt);
##H.tbuf = nan(H.maxlen,1); H.xbuf = nan(H.maxlen, nstate);
##H.head = 1; H.tail = nt; H.count = nt;
##H.tbuf(1:nt) = t_init(:); H.xbuf(1:nt,:) = x_init;
##end
##
##function H = HistoryCache_append(H, t_new, x_new)
##if H.count < H.maxlen
##  H.tail = H.tail + 1; H.tbuf(H.tail) = t_new; H.xbuf(H.tail,:) = x_new(:).'; H.count = H.count + 1;
##else
##  H.head = H.head + 1; if H.head > H.maxlen, H.head = 1; end
##  H.tail = H.head + H.count - 1; if H.tail > H.maxlen, H.tail = H.tail - H.maxlen; end
##  H.tbuf(H.tail) = t_new; H.xbuf(H.tail,:) = x_new(:).';
##end
##end
##
##function xq = HistoryCache_eval(H, tq)
##if H.count == 0, error('HistoryCache is empty'); end
##if H.head <= H.tail
##  tvec = H.tbuf(H.head:H.tail); xmat = H.xbuf(H.head:H.tail,:);
##else
##  tvec = [H.tbuf(H.head:end); H.tbuf(1:H.tail)]; xmat = [H.xbuf(H.head:end,:); H.xbuf(1:H.tail,:)];
##end
##xq = interp1(tvec, xmat, tq, 'pchip')';
##end
##
##%% ------------------------------------------------------------------------
##%% RK4 adaptive (streams accepted steps via callback; returns final state)
##%% ------------------------------------------------------------------------
##function [final_state, stats] = simulate_dde_rk4_adaptive(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
##if nargin < 8, accept_callback = []; end
##at = getfielddef(opts,'atol',1e-6); rt = getfielddef(opts,'rtol',1e-5); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
##interp_grid = getfielddef(opts,'interp_grid',0.05);
##h_eff = max(h_delay,1e-8);
##hist_times = (-h_eff:interp_grid:0);
##hist_states = zeros(numel(hist_times), n_state);
##for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
##H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
##t = tspan(1); x = x0_fun(0); x = x(:);
##stats.rhs_evals = 0; stats.steps = 0;
##h_try = min(max_step, 0.01*(tspan(2)-tspan(1)) + 1e-2);
##while t < tspan(2)-1e-12
##  h_try = min(h_try, tspan(2)-t);
##  x_tau = HistoryCache_eval(H, t - h_delay);
##  k1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
##  x2 = x + 0.5*h_try*k1; x_tau2 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
##  k2 = f_rhs(t + 0.5*h_try, x2, x_tau2); stats.rhs_evals = stats.rhs_evals + 1;
##  x3 = x + 0.5*h_try*k2; x_tau3 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
##  k3 = f_rhs(t + 0.5*h_try, x3, x_tau3); stats.rhs_evals = stats.rhs_evals + 1;
##  x4 = x + h_try*k3; x_tau4 = HistoryCache_eval(H, t + h_try - h_delay);
##  k4 = f_rhs(t + h_try, x4, x_tau4); stats.rhs_evals = stats.rhs_evals + 1;
##  x_rk4 = x + (h_try/6)*(k1 + 2*k2 + 2*k3 + k4);
##  x_heun = x + (h_try/2)*(k1 + k4);
##  err_est = norm(x_rk4 - x_heun, inf);
##  tol = at + rt * max(norm(x,inf), norm(x_rk4,inf));
##  if err_est <= tol
##    t = t + h_try; x = x_rk4; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
##    if ~isempty(accept_callback), accept_callback(t,x); end
##  end
##  if err_est == 0, fac = 2.0; else fac = 0.8 * (tol/err_est)^0.2; end
##  h_try = min(max_step, max(1e-8, fac * h_try));
##end
##final_state = x
##end
##
##%% ------------------------------------------------------------------------
##%% SDIRK(3) (streams accepted steps via callback; returns final state)
##%% ------------------------------------------------------------------------
##function [final_state, stats] = simulate_dde_sdirk(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
##if nargin < 8, accept_callback = []; end
##at = getfielddef(opts,'atol',1e-7); rt = getfielddef(opts,'rtol',1e-7); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
##interp_grid = getfielddef(opts,'interp_grid',0.02); Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
##gamma = 0.435866521508459; A = [gamma 0 0; 0.5-gamma gamma 0; 2*gamma 1-4*gamma gamma]; b = [2*gamma;1-4*gamma;2*gamma];
##h_eff = max(h_delay,1e-8); hist_times = (-h_eff:interp_grid:0); hist_states = zeros(numel(hist_times), n_state);
##for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
##H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
##t = tspan(1); x = x0_fun(0); x = x(:);
##stats.rhs_evals = 0; stats.newton_iters = 0; stats.linear_solves = 0; stats.steps = 0;
##Lmat = extra_args{1}; params = extra_args{2};
##J_prev = []; LU = struct('L',[],'U',[],'P',[],'Q',[],'R',[],'ilu',false);
##h = min(max_step, 1e-2); force_ilu = getfielddef(opts,'force_ilu',false);
##while t < tspan(2)-1e-12
##  h = min(h, tspan(2)-t);
##  K = zeros(n_state,3);
##  for stage = 1:3
##    xi = x;
##    for newt = 1:8
##      ti = t + sum(A(stage,1:stage))*h;
##      x_tau = HistoryCache_eval(H, ti - h_delay);
##      fval = f_rhs(ti, xi, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
##      R = xi - x - h*(A(stage,stage)*fval + sum(A(stage,1:stage-1).*K(:,1:stage-1),2));
##      J = analytic_jacobian_reactor(xi, params, Lmat); M = speye(n_state) - h*A(stage,stage)*J;
##      recompute = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
##      if recompute
##        if force_ilu
##          try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.P=[]; LU.Q=[]; LU.R=[]; LU.ilu=true;
##          catch, [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false; end
##        else
##          try [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false;
##          catch, try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.P=[]; LU.Q=[]; LU.R=[]; LU.ilu=true; catch LU.ilu=false; end
##          end
##        end
##        J_prev = J;
##      end
##      RHS = -R;
##      if LU.ilu
##        try [delta,flag] = gmres(M, RHS, [], 1e-8, 200, LU.L, LU.U); if flag ~= 0, delta = M \ RHS; end
##        catch delta = M \ RHS; end
##      else
##        if ~isempty(LU.R), y = LU.P * RHS; z = LU.L \ y; w = LU.U \ z; w = LU.R \ w; delta = LU.Q * w;
##        else y = LU.P * RHS; z = LU.L \ y; w = LU.U \ z; delta = LU.Q * w; end
##      end
##      xi = xi + delta; stats.newton_iters = stats.newton_iters + 1; stats.linear_solves = stats.linear_solves + 1;
##      if norm(delta, inf) < 1e-8, break; end
##    end
##    K(:,stage) = fval;
##  end
##  x_new = x + h*(K * b);
##  err_est = norm(x_new - x, inf); % conservative
##  tol = at + rt * max(norm(x,inf), norm(x_new,inf));
##  if err_est <= tol
##    t = t + h; x = x_new; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
##    if ~isempty(accept_callback), accept_callback(t,x); end
##    if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
##    h = min(max_step, max(1e-8, fac*h));
##  else
##    h = max(1e-8, 0.5*h);
##  end
##end
##final_state = x
##end
##
##%% ------------------------------------------------------------------------
##%% Rosenbrock (streams accepted steps via callback; returns final state)
##%% ------------------------------------------------------------------------
##function [final_state, stats] = simulate_dde_rosenbrock(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
##if nargin < 8, accept_callback = []; end
##at = getfielddef(opts,'atol',1e-7); rt = getfielddef(opts,'rtol',1e-7); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
##interp_grid = getfielddef(opts,'interp_grid',0.02); Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
##gamma = 0.5;
##h_eff = max(h_delay,1e-8); hist_times = (-h_eff:interp_grid:0); hist_states = zeros(numel(hist_times), n_state);
##for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
##H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
##t = tspan(1); x = x0_fun(0); x = x(:);
##stats.rhs_evals = 0; stats.linear_solves = 0; stats.steps = 0;
##Lmat = extra_args{1}; params = extra_args{2};
##J_prev = []; LU = struct('L',[],'U',[],'P',[],'Q',[],'R',[],'ilu',false);
##h = min(max_step, 1e-2); force_ilu = getfielddef(opts,'force_ilu',false);
##while t < tspan(2)-1e-12
##  h = min(h, tspan(2)-t);
##  x_tau = HistoryCache_eval(H, t - h_delay);
##  J = analytic_jacobian_reactor(x, params, Lmat); M = speye(n_state) - gamma*h*J;
##  recompute = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
##  if recompute
##    if force_ilu
##      try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.ilu=true; catch, [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false; end
##    else
##      try [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false;
##      catch, try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.ilu=true; catch LU.ilu=false; end
##      end
##    end
##    J_prev = J;
##  end
##  f1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
##  if LU.ilu, [k1,fl]=gmres(M,f1,[],1e-8,200,LU.L,LU.U); if fl~=0, k1=M\f1; end
##  else y=LU.P*f1; z=LU.L\y; w=LU.U\z; if ~isempty(LU.R), w=LU.R\w; end; k1=LU.Q*w; end
##  t2 = t + 0.5*h; x2 = x + 0.5*h*k1; x2_tau = HistoryCache_eval(H, t2 - h_delay);
##  f2 = f_rhs(t2, x2, x2_tau); stats.rhs_evals = stats.rhs_evals + 1;
##  rhs2 = f2 + (1/(2*h))*k1;
##  if LU.ilu, [k2,fl]=gmres(M,rhs2,[],1e-8,200,LU.L,LU.U); if fl~=0, k2=M\rhs2; end
##  else y2=LU.P*rhs2; z2=LU.L\y2; w2=LU.U\z2; if ~isempty(LU.R), w2=LU.R\w2; end; k2=LU.Q*w2; end
##  x_new = x + h*(k1 + k2)/2;
##  err_est = norm(x_new - (x + h*k1), inf);
##  tol = at + rt * max(norm(x,inf), norm(x_new,inf));
##  if err_est <= tol
##    t = t + h; x = x_new; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
##    if ~isempty(accept_callback), accept_callback(t,x); end
##    if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
##    h = min(max_step, max(1e-8, fac*h));
##  else
##    h = max(1e-8, 0.5*h);
##  end
##end
##final_state = x
##end
##
##%% ------------------------------------------------------------------------
##%% Plotting from checkpoints (cheap)
##%% ------------------------------------------------------------------------
##function create_3d_plots_from_checkpoints(P, Lmat, Tchk, Xchk, title_prefix)
##if isempty(Tchk) || isempty(Xchk), warning('No checkpoints to plot for %s', title_prefix); return; end
##Nz = P.grid.Nz; Nr = P.grid.Nr; N = P.Ncells;
##cf = Xchk(end,1:N)'; tf = Xchk(end,N+1:end)';
##C = reshape(cf,[Nr Nz]); Tm = reshape(tf,[Nr Nz]);
##z = linspace(0, getfielddef(P,'geometry.L',getfielddef(P,'geometry.tube_length',1.0)), Nz);
##r = linspace(0, getfielddef(P,'geometry.R',getfielddef(P,'geometry.tube_radius',1.0)), Nr);
##max_cells_for_isosurface = 200000;
##if P.Ncells <= max_cells_for_isosurface
##  nang = 24; theta = linspace(0,2*pi,nang);
##  Xcart = zeros(Nr,nang,Nz); Ycart = Xcart; Zcart = zeros(Nr,nang,Nz);
##  Ccart = zeros(Nr,nang,Nz); Tcart = zeros(Nr,nang,Nz);
##  for iz=1:Nz, for ir=1:Nr, for it=1:nang
##    Xcart(ir,it,iz)=r(ir)*cos(theta(it)); Ycart(ir,it,iz)=r(ir)*sin(theta(it)); Zcart(ir,it,iz)=z(iz);
##    Ccart(ir,it,iz)=C(ir,iz); Tcart(ir,it,iz)=Tm(ir,iz);
##  end,end,end
##  figure('Name',[title_prefix ' Temp iso']); p=patch(isosurface(Xcart,Ycart,Zcart,Tcart,mean(Tcart(:)))); isonormals(Xcart,Ycart,Zcart,Tcart,p);
##  set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7); camlight; lighting gouraud; axis equal;
##else
##  figure('Name',[title_prefix ' Conc slice']); imagesc(z,r,C); colorbar; axis xy; xlabel('z (m)'); ylabel('r (m)');
##  figure('Name',[title_prefix ' Centerline T']); plot(z, Tm(ceil(Nr/2),:)); xlabel('z (m)'); ylabel('T (K)');
##end
##end
##
##%% ------------------------------------------------------------------------
##%% Utility
##%% ------------------------------------------------------------------------
##function v = getfielddef(s,f,d)
##if isempty(s) || ~isstruct(s) || ~isfield(s,f), v = d; else v = s.(f); end
##end
##
##
##function industrial_reactors_highres_memory_safe()
##close all force; clc; clearvars -except; rand('state',1);
##
##% ---------------- User geometry and industrial parameters ----------------
##A = set_packed_bed_geometry(12.0,1.6,240,60);
##B = set_multi_tubular_geometry(18.0,0.15,400,32);
##
##A = set_packed_bed_industrial_params(A, struct('Dc',1e-5,'DT',5e-6,'k1',10,'Ea',6e4,'k2',0.5,'beta',5,'gamma',1,'porosity',0.4));
##B = set_multi_tubular_industrial_params(B, struct('Dc',5e-6,'DT',2e-6,'k1',1.2,'Ea',4.5e4,'k2',0.15,'beta',2,'gamma',0.6,'porosity',0.45));
##
##A.Tfinal = 60; B.Tfinal = 120;
##A.solver = default_solver_settings('rosenbrock',1e-7,1e-7,0.02,0.005,1e-8);
##B.solver = default_solver_settings('sdirk',1e-7,1e-7,0.02,0.005,1e-8);
##
##% ---------------- Spatial operators (reuse) ----------------
##L_A = build_axisymmetric_laplacian_kron(A.grid.Nz, A.grid.Nr, A.dx, A.dr)
##L_B = build_axisymmetric_laplacian_kron(B.grid.Nz, B.grid.Nr, B.dx, B.dr)
##
##x0A_fun = @(t) [A.ic.c0; A.ic.T0];
##x0B_fun = @(t) [B.ic.c0; B.ic.T0];
##
##% ---------------- Memory-safe run parameters ----------------
##save_every = 40;            % store every N accepted steps (tune)
##max_nstate_for_lu = 4e4;    % above this prefer ILU+GMRES
##history_margin = 7;        % extra history slots
##
##fprintf('Running high-resolution simulations (memory-safe mode)...\n');
##
##% Reactor A: SDIRK reference (coarse tolerances to reduce cost)
##opts_refA = struct('atol',1e-8,'rtol',1e-8,'max_step',A.Tfinal/500,'interp_grid',A.history.grid,'jacobian_reuse_tol',A.solver.jacobian_reuse_tol);
##[TA_ref_chk, XA_ref_chk, finalA_ref, stats_refA] = run_and_checkpoint(@simulate_dde_sdirk, @(t,x,xt) rhs_packed(t,x,xt,L_A,A), x0A_fun, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {L_A,A}, opts_refA, 'SDIRK-A', save_every, max_nstate_for_lu, history_margin)
##
##% Reactor A: Rosenbrock
##opts_rosA = A.solver.opts; opts_rosA.jacobian_reuse_tol = A.solver.jacobian_reuse_tol;
##[TA_ros_chk, XA_ros_chk, finalA_ros, stats_rosA] = run_and_checkpoint(@simulate_dde_rosenbrock, @(t,x,xt) rhs_packed(t,x,xt,L_A,A), x0A_fun, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {L_A,A}, opts_rosA, 'ROSENBROCK-A', save_every, max_nstate_for_lu, history_margin)
##
##clear TA_ref_chk XA_ref_chk; try pack; catch, end
##
##% Reactor B: SDIRK
##opts_sdirkB = B.solver.opts; opts_sdirkB.jacobian_reuse_tol = B.solver.jacobian_reuse_tol;
##[TB_sdirk_chk, XB_sdirk_chk, finalB_sdirk, stats_sdirkB] = run_and_checkpoint(@simulate_dde_sdirk, @(t,x,xt) rhs_packed(t,x,xt,L_B,B), x0B_fun, B.delay.h, [0 B.Tfinal], 2*B.Ncells, {L_B,B}, opts_sdirkB, 'SDIRK-B', save_every, max_nstate_for_lu, history_margin)
##
##% Save compact results
##save('results_A_rosenbrock.mat','TA_ros_chk','XA_ros_chk','finalA_ros','-v7.3');
##save('results_B_sdirk.mat','TB_sdirk_chk','XB_sdirk_chk','finalB_sdirk','-v7.3');
##
##% Plot from checkpoints (cheap)
##try
##  create_3d_plots_from_checkpoints(A, L_A, TA_ros_chk, XA_ros_chk, 'Packed-bed reactor (A)');
##  create_3d_plots_from_checkpoints(B, L_B, TB_sdirk_chk, XB_sdirk_chk, 'Multi-tubular reactor (B)');
##catch ME
##  warning('Plotting failed: %s', ME.message);
##end
##
##fprintf('All runs complete; checkpoints saved.\n');
##end
##
##
##industrial_reactors_highres_memory_safe()
##

% End of file


#+ more meetings:


##https://payhip.com/b/9Rzy8
##
##https://flippa.com/lps/meeting-scheduler/
##
##https://flippa.com/lps/meeting-scheduler/
##
##So how exactly we will enhance clients' experience on the platform?
##
##Which guidance, which expertise?

% industrial_reactors_highres_memory_safe.m
% Single-file, memory-safe reactor simulator (MATLAB / Octave)
% - Main entry: industrial_reactors_highres_memory_safe()
% - Integrators stream accepted steps to a callback; they do NOT retain full trajectories
% - Checkpoint wrapper collects downsampled checkpoints and final state only
% - Sparse Kronecker Laplacian, analytic Jacobian, LU caching with ILU+GMRES fallback
% Save and run: industrial_reactors_highres_memory_safe()
close all force; clc; clearvars -except; rand('state',1);


%% ------------------------------------------------------------------------
%% run_and_checkpoint: wrapper that collects downsampled checkpoints only
%% ------------------------------------------------------------------------
function [Tchk,Xchk,final_state,stats] = run_and_checkpoint( ...
    solver_fn, f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, ...
    label, save_every, max_nstate_for_lu, history_margin)

% Defaults
if nargin < 11, history_margin = 10; end
if nargin < 10, max_nstate_for_lu = 5e4; end
if nargin < 9, save_every = 50; end

% conservative estimate for preallocation
est_steps = ceil((tspan(2)-tspan(1)) / max(1e-6, getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100))) + 10;
est_store = max(10, ceil(est_steps / max(1, save_every)));

% prepare callback_state in a global (no nested functions)
global run_and_checkpoint_callback_state
run_and_checkpoint_callback_state = struct();
run_and_checkpoint_callback_state.Tchk = zeros(est_store,1);
run_and_checkpoint_callback_state.Xchk = zeros(est_store, n_state);
run_and_checkpoint_callback_state.idx_store = 0;
run_and_checkpoint_callback_state.step_counter = 0;
run_and_checkpoint_callback_state.save_every = save_every;
run_and_checkpoint_callback_state.tspan = tspan;
run_and_checkpoint_callback_state.est_store = est_store;
run_and_checkpoint_callback_state.n_state = n_state;

% adjust history buffer size
interp_grid = getfielddef(opts,'interp_grid',0.01);
opts.history_maxlen = getfielddef(opts,'history_maxlen', ceil(h_delay / interp_grid) + history_margin);

% force ILU for very large systems
if n_state > max_nstate_for_lu
  opts.force_ilu = true;
else
  opts.force_ilu = getfielddef(opts,'force_ilu', false);
end

% Call solver. The solver must accept a callback handle with signature (t_step, x_step)
[final_state, stats] = solver_fn(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, @run_and_checkpoint_callback);

% Extract checkpoints from global state
if isempty(run_and_checkpoint_callback_state) || run_and_checkpoint_callback_state.idx_store == 0
  % No checkpoints were stored by callback: return final state and final time
  Tchk = tspan(2);
  Xchk = final_state(:).';
else
  idx = run_and_checkpoint_callback_state.idx_store;
  Tchk = run_and_checkpoint_callback_state.Tchk(1:idx);
  Xchk = run_and_checkpoint_callback_state.Xchk(1:idx, :);
end

% Clear global callback state to avoid interference with subsequent calls
run_and_checkpoint_callback_state = [];
clear global run_and_checkpoint_callback_state

fprintf('%s done: stored %d checkpoints\n', label, numel(Tchk));
end

% -------------------------------------------------------------------------
% Non-nested callback function used by run_and_checkpoint (must be in same file)
% -------------------------------------------------------------------------
function run_and_checkpoint_callback(t_step, x_step)
global run_and_checkpoint_callback_state
if isempty(run_and_checkpoint_callback_state)
  error('run_and_checkpoint_callback_state is not initialized.');
end

cs = run_and_checkpoint_callback_state;
cs.step_counter = cs.step_counter + 1;

if mod(cs.step_counter, cs.save_every) == 0 || t_step >= cs.tspan(2) - 1e-12
  cs.idx_store = cs.idx_store + 1;
  idx = cs.idx_store;
  if idx > size(cs.Tchk,1)
    extra = cs.est_store;
    cs.Tchk = [cs.Tchk; zeros(extra,1)];
    cs.Xchk = [cs.Xchk; zeros(extra, cs.n_state)];
  end
  cs.Tchk(idx) = t_step;
  cs.Xchk(idx, :) = x_step(:).';
end

% write back to global
run_and_checkpoint_callback_state = cs;
end


%% ------------------------------------------------------------------------
%% Geometry setters
%% ------------------------------------------------------------------------
function R = set_packed_bed_geometry(L, Rrad, Nz, Nr)
R = struct();
R.geometry.L = L; R.geometry.R = Rrad;
R.grid.Nz = Nz; R.grid.Nr = Nr; R.Ncells = Nz * Nr;
R.dx = L / (Nz - 1); R.dr = Rrad / (Nr - 1);
R.ic.c0 = 0.5 * ones(R.Ncells,1); R.ic.T0 = 300.0 * ones(R.Ncells,1);
R.delay.h = 0.3; R.history.grid = 0.005;
R.bc.inlet.c = 0.8; R.bc.inlet.T = 300.0;
end

function R = set_multi_tubular_geometry(Ltube, Rtube, Nz, Nr)
R = struct();
R.geometry.tube_length = Ltube; R.geometry.tube_radius = Rtube;
R.grid.Nz = Nz; R.grid.Nr = Nr; R.Ncells = Nz * Nr;
R.dx = Ltube / (Nz - 1); R.dr = Rtube / (Nr - 1);
R.ic.c0 = 0.4 * ones(R.Ncells,1); R.ic.T0 = 320.0 * ones(R.Ncells,1);
R.delay.h = 0.15; R.history.grid = 0.005;
R.bc.inlet.c = 0.6; R.bc.inlet.T = 320.0;
end

%% ------------------------------------------------------------------------
%% Industrial params with defaults
%% ------------------------------------------------------------------------
function s = default_solver_settings(method, atol, rtol, max_step, interp_grid, jacobian_reuse_tol)
if nargin < 1, method = 'rosenbrock'; end
if nargin < 2, atol = 1e-7; end
if nargin < 3, rtol = 1e-7; end
if nargin < 4, max_step = 0.02; end
if nargin < 5, interp_grid = 0.005; end
if nargin < 6, jacobian_reuse_tol = 1e-8; end
s = struct(); s.method = method; s.opts = struct('atol',atol,'rtol',rtol,'max_step',max_step,'interp_grid',interp_grid); s.jacobian_reuse_tol = jacobian_reuse_tol;
end

function R = set_packed_bed_industrial_params(R, phys)
if nargin < 2, phys = struct(); end
phys.Dc = getfielddef(phys,'Dc',1.0e-5);
phys.DT = getfielddef(phys,'DT',5.0e-6);
phys.k1 = getfielddef(phys,'k1',10.0);
phys.Ea = getfielddef(phys,'Ea',6.0e4);
phys.k2 = getfielddef(phys,'k2',0.5);
phys.beta = getfielddef(phys,'beta',5.0);
phys.gamma = getfielddef(phys,'gamma',1.0);
phys.porosity = getfielddef(phys,'porosity',0.4);
phys.rho = getfielddef(phys,'rho',1.2);
phys.cp = getfielddef(phys,'cp',1000.0);
R.phys = phys;
end

function R = set_multi_tubular_industrial_params(R, phys)
if nargin < 2, phys = struct(); end
phys.Dc = getfielddef(phys,'Dc',5.0e-6);
phys.DT = getfielddef(phys,'DT',2.0e-6);
phys.k1 = getfielddef(phys,'k1',1.2);
phys.Ea = getfielddef(phys,'Ea',4.5e4);
phys.k2 = getfielddef(phys,'k2',0.15);
phys.beta = getfielddef(phys,'beta',2.0);
phys.gamma = getfielddef(phys,'gamma',0.6);
phys.porosity = getfielddef(phys,'porosity',0.45);
phys.rho = getfielddef(phys,'rho',1.0);
phys.cp = getfielddef(phys,'cp',1200.0);
R.phys = phys;
end

%% ------------------------------------------------------------------------
%% Axisymmetric Laplacian (Kronecker)
%% ------------------------------------------------------------------------
function L = build_axisymmetric_laplacian_kron(Nz, Nr, dz, dr)
ex = spdiags([ones(Nz,1), -2*ones(Nz,1), ones(Nz,1)], [-1,0,1], Nz, Nz);
ey = spdiags([ones(Nr,1), -2*ones(Nr,1), ones(Nr,1)], [-1,0,1], Nr, Nr);
Ix = speye(Nz); Iy = speye(Nr);
Lz = ex / (dz*dz); Lr = ey / (dr*dr);
L = kron(Iy, Lz) + kron(Lr, Ix);
end

%% ------------------------------------------------------------------------
%% Reactor RHS (packed)
%% ------------------------------------------------------------------------
function dX = rhs_packed(t, X, X_tau, Lmat, P)
N = P.Ncells;
c = X(1:N); T = X(N+1:end);
Rrate = P.phys.k1 .* c .* exp(-(P.phys.Ea/8.314) .* (1./T - 1/300));
dc = P.phys.Dc * (Lmat * c) - Rrate - P.phys.k2 .* c;
dT = P.phys.DT * (Lmat * T) + P.phys.beta .* Rrate - P.phys.gamma .* (X_tau(N+1:end) - P.bc.inlet.T);
dX = [dc; dT];
end

%% ------------------------------------------------------------------------
%% Analytic Jacobian (sparse block)
%% ------------------------------------------------------------------------
function J = analytic_jacobian_reactor(x, params, Lmat)
N = params.Ncells;
c = x(1:N); T = x(N+1:end);
expT = exp(-(params.phys.Ea/8.314) .* (1./T - 1/300));
dr_dc = params.phys.k1 .* expT;
dexp_dT = expT .* (params.phys.Ea/8.314) .* (1./(T.^2));
dr_dT = params.phys.k1 .* c .* dexp_dT;
J11 = params.phys.Dc * Lmat - spdiags(params.phys.k2 * ones(N,1) + dr_dc, 0, N, N);
J12 = - spdiags(dr_dT, 0, N, N);
J21 = spdiags(params.phys.beta * dr_dc, 0, N, N);
J22 = params.phys.DT * Lmat - spdiags(params.phys.gamma * ones(N,1), 0, N, N);
J = [J11, J12; J21, J22];
end

%% ------------------------------------------------------------------------
%% HistoryCache (circular)
%% ------------------------------------------------------------------------
function H = HistoryCache_create(t_init, x_init, maxlen)
if nargin < 3, maxlen = 200000; end
nt = numel(t_init); nstate = size(x_init,2);
H.maxlen = max(maxlen, nt);
H.tbuf = nan(H.maxlen,1); H.xbuf = nan(H.maxlen, nstate);
H.head = 1; H.tail = nt; H.count = nt;
H.tbuf(1:nt) = t_init(:); H.xbuf(1:nt,:) = x_init;
end

function H = HistoryCache_append(H, t_new, x_new)
if H.count < H.maxlen
  H.tail = H.tail + 1; H.tbuf(H.tail) = t_new; H.xbuf(H.tail,:) = x_new(:).'; H.count = H.count + 1;
else
  H.head = H.head + 1; if H.head > H.maxlen, H.head = 1; end
  H.tail = H.head + H.count - 1; if H.tail > H.maxlen, H.tail = H.tail - H.maxlen; end
  H.tbuf(H.tail) = t_new; H.xbuf(H.tail,:) = x_new(:).';
end
end

function xq = HistoryCache_eval(H, tq)
if H.count == 0, error('HistoryCache is empty'); end
if H.head <= H.tail
  tvec = H.tbuf(H.head:H.tail); xmat = H.xbuf(H.head:H.tail,:);
else
  tvec = [H.tbuf(H.head:end); H.tbuf(1:H.tail)]; xmat = [H.xbuf(H.head:end,:); H.xbuf(1:H.tail,:)];
end
xq = interp1(tvec, xmat, tq, 'pchip')';
end

%% ------------------------------------------------------------------------
%% RK4 adaptive (streams accepted steps via callback; returns final state)
%% ------------------------------------------------------------------------
function [final_state, stats] = simulate_dde_rk4_adaptive(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
if nargin < 8, accept_callback = []; end
at = getfielddef(opts,'atol',1e-6); rt = getfielddef(opts,'rtol',1e-5); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
interp_grid = getfielddef(opts,'interp_grid',0.05);
h_eff = max(h_delay,1e-8);
hist_times = (-h_eff:interp_grid:0);
hist_states = zeros(numel(hist_times), n_state);
for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
t = tspan(1); x = x0_fun(0); x = x(:);
stats.rhs_evals = 0; stats.steps = 0;
h_try = min(max_step, 0.01*(tspan(2)-tspan(1)) + 1e-2);
while t < tspan(2)-1e-12
  h_try = min(h_try, tspan(2)-t);
  x_tau = HistoryCache_eval(H, t - h_delay);
  k1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
  x2 = x + 0.5*h_try*k1; x_tau2 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
  k2 = f_rhs(t + 0.5*h_try, x2, x_tau2); stats.rhs_evals = stats.rhs_evals + 1;
  x3 = x + 0.5*h_try*k2; x_tau3 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
  k3 = f_rhs(t + 0.5*h_try, x3, x_tau3); stats.rhs_evals = stats.rhs_evals + 1;
  x4 = x + h_try*k3; x_tau4 = HistoryCache_eval(H, t + h_try - h_delay);
  k4 = f_rhs(t + h_try, x4, x_tau4); stats.rhs_evals = stats.rhs_evals + 1;
  x_rk4 = x + (h_try/6)*(k1 + 2*k2 + 2*k3 + k4);
  x_heun = x + (h_try/2)*(k1 + k4);
  err_est = norm(x_rk4 - x_heun, inf);
  tol = at + rt * max(norm(x,inf), norm(x_rk4,inf));
  if err_est <= tol
    t = t + h_try; x = x_rk4; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
    if ~isempty(accept_callback), accept_callback(t,x); end
  end
  if err_est == 0, fac = 2.0; else fac = 0.8 * (tol/err_est)^0.2; end
  h_try = min(max_step, max(1e-8, fac * h_try));
end
final_state = x;
end

%% ------------------------------------------------------------------------
%% SDIRK(3) (streams accepted steps via callback; returns final state)
%% ------------------------------------------------------------------------
function [final_state, stats] = simulate_dde_sdirk(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
if nargin < 8, accept_callback = []; end
at = getfielddef(opts,'atol',1e-7); rt = getfielddef(opts,'rtol',1e-7); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
interp_grid = getfielddef(opts,'interp_grid',0.02); Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
gamma = 0.435866521508459; A = [gamma 0 0; 0.5-gamma gamma 0; 2*gamma 1-4*gamma gamma]; b = [2*gamma;1-4*gamma;2*gamma];
h_eff = max(h_delay,1e-8); hist_times = (-h_eff:interp_grid:0); hist_states = zeros(numel(hist_times), n_state);
for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
t = tspan(1); x = x0_fun(0); x = x(:);
stats.rhs_evals = 0; stats.newton_iters = 0; stats.linear_solves = 0; stats.steps = 0;
Lmat = extra_args{1}; params = extra_args{2};
J_prev = []; LU = struct('L',[],'U',[],'P',[],'Q',[],'R',[],'ilu',false);
h = min(max_step, 1e-2); force_ilu = getfielddef(opts,'force_ilu',false);
while t < tspan(2)-1e-12
  h = min(h, tspan(2)-t); K = zeros(n_state,3);
  for stage = 1:3
    xi = x;
    for newt = 1:8
      ti = t + sum(A(stage,1:stage))*h;
      x_tau = HistoryCache_eval(H, ti - h_delay);
      fval = f_rhs(ti, xi, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
      R = xi - x - h*(A(stage,stage)*fval + sum(A(stage,1:stage-1).*K(:,1:stage-1),2));
      J = analytic_jacobian_reactor(xi, params, Lmat); M = speye(n_state) - h*A(stage,stage)*J;
      recompute = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
      if recompute
        if force_ilu
          try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.P=[]; LU.Q=[]; LU.R=[]; LU.ilu=true;
          catch, [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false; end
        else
          try [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false;
          catch, try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.P=[]; LU.Q=[]; LU.R=[]; LU.ilu=true; catch LU.ilu=false; end
          end
        end
        J_prev = J;
      end
      RHS = -R;
      if LU.ilu
        try [delta,flag] = gmres(M, RHS, [], 1e-8, 200, LU.L, LU.U); if flag ~= 0, delta = M \ RHS; end
        catch delta = M \ RHS; end
      else
        if ~isempty(LU.R), y = LU.P * RHS; z = LU.L \ y; w = LU.U \ z; w = LU.R \ w; delta = LU.Q * w;
        else y = LU.P * RHS; z = LU.L \ y; w = LU.U \ z; delta = LU.Q * w; end
      end
      xi = xi + delta; stats.newton_iters = stats.newton_iters + 1; stats.linear_solves = stats.linear_solves + 1;
      if norm(delta, inf) < 1e-8, break; end
    end
    K(:,stage) = fval;
  end
  x_new = x + h*(K * b);
  err_est = norm(x_new - x, inf); % conservative
  tol = at + rt * max(norm(x,inf), norm(x_new,inf));
  if err_est <= tol
    t = t + h; x = x_new; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
    if ~isempty(accept_callback), accept_callback(t,x); end
    if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
    h = min(max_step, max(1e-8, fac*h));
  else
    h = max(1e-8, 0.5*h);
  end
end
final_state = x;
end

%% ------------------------------------------------------------------------
%% Rosenbrock (streams accepted steps via callback; returns final state)
%% ------------------------------------------------------------------------
function [final_state, stats] = simulate_dde_rosenbrock(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, accept_callback)
if nargin < 8, accept_callback = []; end
at = getfielddef(opts,'atol',1e-7); rt = getfielddef(opts,'rtol',1e-7); max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
interp_grid = getfielddef(opts,'interp_grid',0.02); Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
gamma = 0.5;
h_eff = max(h_delay,1e-8); hist_times = (-h_eff:interp_grid:0); hist_states = zeros(numel(hist_times), n_state);
for k=1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
H = HistoryCache_create(hist_times, hist_states, getfielddef(opts,'history_maxlen',200000));
t = tspan(1); x = x0_fun(0); x = x(:);
stats.rhs_evals = 0; stats.linear_solves = 0; stats.steps = 0;
Lmat = extra_args{1}; params = extra_args{2};
J_prev = []; LU = struct('L',[],'U',[],'P',[],'Q',[],'R',[],'ilu',false);
h = min(max_step, 1e-2); force_ilu = getfielddef(opts,'force_ilu',false);
while t < tspan(2)-1e-12
  h = min(h, tspan(2)-t);
  x_tau = HistoryCache_eval(H, t - h_delay);
  J = analytic_jacobian_reactor(x, params, Lmat); M = speye(n_state) - gamma*h*J;
  recompute = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
  if recompute
    if force_ilu
      try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.ilu=true; catch, [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false; end
    else
      try [Lfac,Ufac,Pfac,Qfac,Rfac]=lu(M); LU.L=Lfac; LU.U=Ufac; LU.P=Pfac; LU.Q=Qfac; LU.R=Rfac; LU.ilu=false;
      catch, try st.type='ilutp'; st.droptol=1e-3; [Lilu,Uilu]=ilu(M,st); LU.L=Lilu; LU.U=Uilu; LU.ilu=true; catch LU.ilu=false; end
      end
    end
    J_prev = J;
  end
  f1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
  if LU.ilu, [k1,fl]=gmres(M,f1,[],1e-8,200,LU.L,LU.U); if fl~=0, k1=M\f1; end
  else y=LU.P*f1; z=LU.L\y; w=LU.U\z; if ~isempty(LU.R), w=LU.R\w; end; k1=LU.Q*w; end
  t2 = t + 0.5*h; x2 = x + 0.5*h*k1; x2_tau = HistoryCache_eval(H, t2 - h_delay);
  f2 = f_rhs(t2, x2, x2_tau); stats.rhs_evals = stats.rhs_evals + 1;
  rhs2 = f2 + (1/(2*h))*k1;
  if LU.ilu, [k2,fl]=gmres(M,rhs2,[],1e-8,200,LU.L,LU.U); if fl~=0, k2=M\rhs2; end
  else y2=LU.P*rhs2; z2=LU.L\y2; w2=LU.U\z2; if ~isempty(LU.R), w2=LU.R\w2; end; k2=LU.Q*w2; end
  x_new = x + h*(k1 + k2)/2;
  err_est = norm(x_new - (x + h*k1), inf);
  tol = at + rt * max(norm(x,inf), norm(x_new,inf));
  if err_est <= tol
    t = t + h; x = x_new; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
    if ~isempty(accept_callback), accept_callback(t,x); end
    if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
    h = min(max_step, max(1e-8, fac*h));
  else
    h = max(1e-8, 0.5*h);
  end
end
final_state = x;
end

%% ------------------------------------------------------------------------
%% Plotting from checkpoints (cheap)
%% ------------------------------------------------------------------------
function create_3d_plots_from_checkpoints(P, Lmat, Tchk, Xchk, title_prefix)
if isempty(Tchk) || isempty(Xchk), warning('No checkpoints to plot for %s', title_prefix); return; end
Nz = P.grid.Nz; Nr = P.grid.Nr; N = P.Ncells;
cf = Xchk(end,1:N)'; tf = Xchk(end,N+1:end)';
C = reshape(cf,[Nr Nz]); Tm = reshape(tf,[Nr Nz]);
z = linspace(0, getfielddef(P,'geometry.L',getfielddef(P,'geometry.tube_length',1.0)), Nz);
r = linspace(0, getfielddef(P,'geometry.R',getfielddef(P,'geometry.tube_radius',1.0)), Nr);
max_cells_for_isosurface = 200000;
if P.Ncells <= max_cells_for_isosurface
  nang = 24; theta = linspace(0,2*pi,nang);
  Xcart = zeros(Nr,nang,Nz); Ycart = Xcart; Zcart = zeros(Nr,nang,Nz);
  Ccart = zeros(Nr,nang,Nz); Tcart = zeros(Nr,nang,Nz);
  for iz=1:Nz, for ir=1:Nr, for it=1:nang
    Xcart(ir,it,iz)=r(ir)*cos(theta(it)); Ycart(ir,it,iz)=r(ir)*sin(theta(it)); Zcart(ir,it,iz)=z(iz);
    Ccart(ir,it,iz)=C(ir,iz); Tcart(ir,it,iz)=Tm(ir,iz);
  end,end,end
  figure('Name',[title_prefix ' Temp iso']); p=patch(isosurface(Xcart,Ycart,Zcart,Tcart,mean(Tcart(:)))); isonormals(Xcart,Ycart,Zcart,Tcart,p);
  set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7); camlight; lighting gouraud; axis equal;
else
  figure('Name',[title_prefix ' Conc slice']); imagesc(z,r,C); colorbar; axis xy; xlabel('z (m)'); ylabel('r (m)');
  figure('Name',[title_prefix ' Centerline T']); plot(z, Tm(ceil(Nr/2),:)); xlabel('z (m)'); ylabel('T (K)');
end
end

%% ------------------------------------------------------------------------
%% Utility
%% ------------------------------------------------------------------------
function v = getfielddef(s,f,d)
if isempty(s) || ~isstruct(s) || ~isfield(s,f), v = d; else v = s.(f); end
end



function industrial_reactors_highres_memory_safe()

% ---------------- User geometry and industrial parameters ----------------
A = set_packed_bed_geometry(12.0,1.6,240,60);
B = set_multi_tubular_geometry(18.0,0.15,400,32);

A = set_packed_bed_industrial_params(A, struct('Dc',1e-5,'DT',5e-6,'k1',10,'Ea',6e4,'k2',0.5,'beta',5,'gamma',1,'porosity',0.4));
B = set_multi_tubular_industrial_params(B, struct('Dc',5e-6,'DT',2e-6,'k1',1.2,'Ea',4.5e4,'k2',0.15,'beta',2,'gamma',0.6,'porosity',0.45));

A.Tfinal = 60; B.Tfinal = 120;
A.solver = default_solver_settings('rosenbrock',1e-7,1e-7,0.02,0.005,1e-8);
B.solver = default_solver_settings('sdirk',1e-7,1e-7,0.02,0.005,1e-8);

% ---------------- Spatial operators (reuse) ----------------
L_A = build_axisymmetric_laplacian_kron(A.grid.Nz, A.grid.Nr, A.dx, A.dr);
L_B = build_axisymmetric_laplacian_kron(B.grid.Nz, B.grid.Nr, B.dx, B.dr);

x0A_fun = @(t) [A.ic.c0; A.ic.T0];
x0B_fun = @(t) [B.ic.c0; B.ic.T0];

% ---------------- Memory-safe run parameters ----------------
save_every = 50;            % store every N accepted steps (tune)
max_nstate_for_lu = 5e4;    % above this prefer ILU+GMRES
history_margin = 10;        % extra history slots

fprintf('Running high-resolution simulations (memory-safe mode)...\n');

% Reactor A: SDIRK reference (coarse tolerances to reduce cost)
opts_refA = struct('atol',1e-8,'rtol',1e-8,'max_step',A.Tfinal/500,'interp_grid',A.history.grid,'jacobian_reuse_tol',A.solver.jacobian_reuse_tol);
[TA_ref_chk, XA_ref_chk, finalA_ref, stats_refA] = run_and_checkpoint(@simulate_dde_sdirk, @(t,x,xt) rhs_packed(t,x,xt,L_A,A), x0A_fun, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {L_A,A}, opts_refA, 'SDIRK-A', save_every, max_nstate_for_lu, history_margin);

% Reactor A: Rosenbrock
opts_rosA = A.solver.opts; opts_rosA.jacobian_reuse_tol = A.solver.jacobian_reuse_tol;
[TA_ros_chk, XA_ros_chk, finalA_ros, stats_rosA] = run_and_checkpoint(@simulate_dde_rosenbrock, @(t,x,xt) rhs_packed(t,x,xt,L_A,A), x0A_fun, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {L_A,A}, opts_rosA, 'ROSENBROCK-A', save_every, max_nstate_for_lu, history_margin);

clear TA_ref_chk XA_ref_chk; try pack; catch, end

% Reactor B: SDIRK
opts_sdirkB = B.solver.opts; opts_sdirkB.jacobian_reuse_tol = B.solver.jacobian_reuse_tol;
[TB_sdirk_chk, XB_sdirk_chk, finalB_sdirk, stats_sdirkB] = run_and_checkpoint(@simulate_dde_sdirk, @(t,x,xt) rhs_packed(t,x,xt,L_B,B), x0B_fun, B.delay.h, [0 B.Tfinal], 2*B.Ncells, {L_B,B}, opts_sdirkB, 'SDIRK-B', save_every, max_nstate_for_lu, history_margin);

% Save compact results
save('results_A_rosenbrock.mat','TA_ros_chk','XA_ros_chk','finalA_ros','-v7.3');
save('results_B_sdirk.mat','TB_sdirk_chk','XB_sdirk_chk','finalB_sdirk','-v7.3');

% Plot from checkpoints (cheap)
try
  create_3d_plots_from_checkpoints(A, L_A, TA_ros_chk, XA_ros_chk, 'Packed-bed reactor (A)');
  create_3d_plots_from_checkpoints(B, L_B, TB_sdirk_chk, XB_sdirk_chk, 'Multi-tubular reactor (B)');
catch ME
  warning('Plotting failed: %s', ME.message);
end

fprintf('All runs complete; checkpoints saved.\n');
end

% End of file

industrial_reactors_highres_memory_safe()



###Prod. by Dipl.-Ing. Mr. Paul Mazniker + MSCopilot
###  +380990535261   https://diag.net/u/u6r3ondjie0w0l8138bafm095b
###   https://github.com/goodengineer
###   Private Zoom Meeting: https://us05web.zoom.us/j/7345311127?pwd=k7GWLjDGkVM3Ntpa6UqIXsBD8BPbGZ.1
###   Private GMeet Meeting: https://meet.google.com/ays-zskm-yxv
###   https://gust.com/companies/mr-pavel-s-startup
###   https://orcid.org/0000-0001-8184-8166
###   https://willwork781147312.wordpress.com/portfolio/cp/
###   https://www.mathworks.com/matlabcentral/profile/authors/1643504
###   https://www.researchgate.net/profile/Pavel-Mazniker
###   https://nanohub.org/members/130066
###   ko-fi.com/mrpaul61277
###   https://payhip.com/b/9Rzy8
###   https://independant.academia.edu/PavelMazniker
###   https://www.youtube.com/channel/UCC__7jMOAHak0MVkUFtmO-w
###   https://scholar.google.co.uk/citations?user=cnTX270AAAAJ&hl=en
###   https://flippa.com/users/4752309
###   https://wellfound.com/u/mr-p-mazniker


