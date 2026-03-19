% industrial_reactors_highres_full_fixed.m
% Consolidated, fixed, and optimized version of industrial_reactors_highres_full.m
% - Robust sparse LU usage with permutation-aware solves and ILU+GMRES fallback
% - HistoryCache implemented as a circular buffer (preallocated)
% - Analytic Jacobian built as sparse blocks (no dense intermediates)
% - LU factorization caching and reuse with tolerance
% - Plotting optimized to avoid expensive isosurface on huge grids
% - Minor vectorization and preallocation improvements
% Usage:
%   Add to MATLAB/Octave path and run:
%     industrial_reactors_highres_full_fixed
%
% Note: For Octave, ensure required packages (sparse) are loaded.
  close all force;
  clc;
  clear all
  clear global;
  clear classes;

  rand('state',1);

  % Environment tuning (Octave-friendly guards)
  try
    pack;
  catch
  end

  %% ------------------------------------------------------------------------
%% Geometry and parameter setters (unchanged)
%% ------------------------------------------------------------------------
function R = set_packed_bed_geometry(L, Rrad, Nz, Nr)
  R = struct();
  R.geometry.L = L;
  R.geometry.R = Rrad;
  R.grid.Nz = Nz;
  R.grid.Nr = Nr;
  R.Ncells = Nz * Nr;
  R.dx = L / (Nz - 1);
  R.dr = Rrad / (Nr - 1);
  R.ic.c0 = 0.5 * ones(R.Ncells,1);
  R.ic.T0 = 300.0 * ones(R.Ncells,1);
  R.delay.h = 0.3;
  R.history.grid = 0.005;
  R.bc.inlet.c = 0.8;
  R.bc.inlet.T = 300.0;
end

function R = set_multi_tubular_geometry(Ltube, Rtube, Nz, Nr)
  R = struct();
  R.geometry.tube_length = Ltube;
  R.geometry.tube_radius = Rtube;
  R.grid.Nz = Nz;
  R.grid.Nr = Nr;
  R.Ncells = Nz * Nr;
  R.dx = Ltube / (Nz - 1);
  R.dr = Rtube / (Nr - 1);
  R.ic.c0 = 0.4 * ones(R.Ncells,1);
  R.ic.T0 = 320.0 * ones(R.Ncells,1);
  R.delay.h = 0.15;
  R.history.grid = 0.005;
  R.bc.inlet.c = 0.6;
  R.bc.inlet.T = 320.0;
end

function R = set_packed_bed_industrial_params(R, phys)
  R.phys = phys;
  R.phys.rho = getfielddef(phys,'rho',1.2);
  R.phys.cp  = getfielddef(phys,'cp',1000);
end

function R = set_multi_tubular_industrial_params(R, phys)
  R.phys = phys;
  R.phys.rho = getfielddef(phys,'rho',1.0);
  R.phys.cp  = getfielddef(phys,'cp',1200);
end

function s = default_solver_settings(method, atol, rtol, max_step, interp_grid, jacobian_reuse_tol)
  s = struct();
  s.method = method;
  s.opts = struct('atol',atol,'rtol',rtol,'max_step',max_step,'interp_grid',interp_grid);
  s.jacobian_reuse_tol = jacobian_reuse_tol;
end

%% ------------------------------------------------------------------------
%% Spatial operator: axisymmetric Laplacian (z-r) -> 1D indexing
%% ------------------------------------------------------------------------
function L = build_axisymmetric_laplacian(Nz, Nr, dz, dr)
  N = Nz * Nr;
  L = spalloc(N,N,5*N);
  idx = @(iz,ir) (ir-1)*Nz + iz;
  for ir = 1:Nr
    for iz = 1:Nz
      k = idx(iz,ir);
      if iz > 1, L(k, idx(iz-1,ir)) = 1/dz^2; end
      if iz < Nz, L(k, idx(iz+1,ir)) = 1/dz^2; end
      if ir > 1, L(k, idx(iz,ir-1)) = 1/dr^2; end
      if ir < Nr, L(k, idx(iz,ir+1)) = 1/dr^2; end
      L(k,k) = -sum(L(k,:));
    end
  end
end

%% ------------------------------------------------------------------------
%% Reactor RHS (shared form). Expects X = [c; T], X_tau same shape.
%% ------------------------------------------------------------------------
function dX = rhs_packed(t, X, X_tau, Lmat, P)
  N = P.Ncells;
  c = X(1:N); T = X(N+1:end);
  R = P.phys.k1 .* c .* exp(-(P.phys.Ea/8.314) .* (1./T - 1/300));
  dc = P.phys.Dc * (Lmat * c) - R - P.phys.k2 .* c;
  dT = P.phys.DT * (Lmat * T) + P.phys.beta .* R - P.phys.gamma .* (X_tau(N+1:end) - P.bc.inlet.T);
  dX = [dc; dT];
end

%% ------------------------------------------------------------------------
%% History cache (circular buffer implementation)
%% ------------------------------------------------------------------------
function H = HistoryCache_create(t_init, x_init, maxlen)
  if nargin < 3, maxlen = 200000; end
  nt = numel(t_init);
  nstate = size(x_init,2);
  H.maxlen = max(maxlen, nt);
  H.tbuf = nan(H.maxlen,1);
  H.xbuf = nan(H.maxlen, nstate);
  H.head = 1; H.tail = nt; H.count = nt;
  H.tbuf(1:nt) = t_init(:);
  H.xbuf(1:nt,:) = x_init;
end

function H = HistoryCache_append(H, t_new, x_new)
  if H.count < H.maxlen
    H.tail = H.tail + 1;
    H.tbuf(H.tail) = t_new;
    H.xbuf(H.tail,:) = x_new(:).';
    H.count = H.count + 1;
  else
    H.head = H.head + 1;
    if H.head > H.maxlen, H.head = 1; end
    H.tail = H.head + H.count - 1;
    if H.tail > H.maxlen, H.tail = H.tail - H.maxlen; end
    H.tbuf(H.tail) = t_new;
    H.xbuf(H.tail,:) = x_new(:).';
  end
end

function xq = HistoryCache_eval(H, tq)
  if H.count == 0, error('HistoryCache is empty'); end
  if H.head <= H.tail
    tvec = H.tbuf(H.head:H.tail);
    xmat = H.xbuf(H.head:H.tail,:);
  else
    tvec = [H.tbuf(H.head:end); H.tbuf(1:H.tail)];
    xmat = [H.xbuf(H.head:end,:); H.xbuf(1:H.tail,:)];
  end
  if isscalar(tq)
    xq = interp1(tvec, xmat, tq, 'pchip')';
  else
    xq = interp1(tvec, xmat, tq, 'pchip');
  end
end

%% ------------------------------------------------------------------------
%% Analytic Jacobian builder (sparse block assembly)
%% ------------------------------------------------------------------------
function J = analytic_jacobian_reactor(x, params, Lmat)
  N = params.Ncells;
  c = x(1:N); T = x(N+1:end);
  k1 = params.phys.k1; k2 = params.phys.k2; beta = params.phys.beta; gamma = params.phys.gamma;
  expT = exp(-(params.phys.Ea/8.314) .* (1./T - 1/300));
  diag_dr_dc = k1 .* expT;
  dexp_dT = expT .* (params.phys.Ea/8.314) .* (1./(T.^2));
  diag_dr_dT = k1 .* c .* dexp_dT;
  J11 = params.phys.Dc * Lmat - spdiags(k2 * ones(N,1) + diag_dr_dc, 0, N, N);
  J12 = - spdiags(diag_dr_dT, 0, N, N);
  J21 = spdiags(beta * diag_dr_dc, 0, N, N);
  J22 = params.phys.DT * Lmat - spdiags(gamma * ones(N,1), 0, N, N);
  J = [J11, J12; J21, J22];
end

%% ------------------------------------------------------------------------
%% Integrators (RK4 adaptive, SDIRK with robust LU caching, Rosenbrock)
%% ------------------------------------------------------------------------
function [T_out, X_out, stats] = simulate_dde_rk4_adaptive(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts)
  if nargin < 7, opts = struct(); end
  atol = getfielddef(opts,'atol',1e-4);
  rtol = getfielddef(opts,'rtol',1e-3);
  max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
  interp_grid = getfielddef(opts,'interp_grid',0.05);

  h_eff = max(h_delay, 1e-8);
  hist_times = (-h_eff:interp_grid:0);
  hist_states = zeros(numel(hist_times), n_state);
  for k = 1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
  H = HistoryCache_create(hist_times, hist_states, 200000);

  t = tspan(1); x = x0_fun(0); x = x(:);
  T_out = t; X_out = x.'; stats.rhs_evals = 0; stats.steps = 0;
  h_try = min(max_step, 0.01*(tspan(2)-tspan(1)) + 1e-2);

  while t < tspan(2) - 1e-12
    h_try = min(h_try, tspan(2) - t);
    x_tau = HistoryCache_eval(H, t - h_delay);
    k1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
    x2 = x + 0.5*h_try*k1;
    x_tau2 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
    k2 = f_rhs(t + 0.5*h_try, x2, x_tau2); stats.rhs_evals = stats.rhs_evals + 1;
    x3 = x + 0.5*h_try*k2;
    x_tau3 = HistoryCache_eval(H, t + 0.5*h_try - h_delay);
    k3 = f_rhs(t + 0.5*h_try, x3, x_tau3); stats.rhs_evals = stats.rhs_evals + 1;
    x4 = x + h_try*k3;
    x_tau4 = HistoryCache_eval(H, t + h_try - h_delay);
    k4 = f_rhs(t + h_try, x4, x_tau4); stats.rhs_evals = stats.rhs_evals + 1;

    x_rk4 = x + (h_try/6)*(k1 + 2*k2 + 2*k3 + k4);
    x_heun = x + (h_try/2)*(k1 + k4);
    err_est = norm(x_rk4 - x_heun, inf);
    tol = atol + rtol * max(norm(x,inf), norm(x_rk4,inf));

    if err_est <= tol
      t = t + h_try; x = x_rk4;
      T_out(end+1,1) = t; X_out(end+1,:) = x.'; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
    end
    if err_est == 0, fac = 2.0; else fac = 0.8 * (tol/err_est)^0.2; end
    h_try = min(max_step, max(1e-6, fac * h_try))
  end
end

function [T_out, X_out, stats] = simulate_dde_sdirk(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts)
  if nargin < 7, opts = struct(); end
  atol = getfielddef(opts,'atol',1e-5);
  rtol = getfielddef(opts,'rtol',1e-5);
  max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
  interp_grid = getfielddef(opts,'interp_grid',0.02);

  gamma = 0.435866521508458999416019;
  A = [gamma, 0, 0; 0.5-gamma, gamma, 0; 2*gamma, 1-4*gamma, gamma];
  b = [2*gamma; 1-4*gamma; 2*gamma];

  h_eff = max(h_delay, 1e-8);
  hist_times = (-h_eff:interp_grid:0);
  hist_states = zeros(numel(hist_times), n_state);
  for k = 1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
  H = HistoryCache_create(hist_times, hist_states, 200000);

  t = tspan(1); x = x0_fun(0); x = x(:);
  T_out = t; X_out = x.'; stats.rhs_evals = 0; stats.newton_iters = 0; stats.linear_solves = 0; stats.steps = 0;

  Lmat = []; params = [];
  if ~isempty(extra_args)
    for k = 1:numel(extra_args)
      if issparse(extra_args{k}) && isempty(Lmat), Lmat = extra_args{k};
      elseif isstruct(extra_args{k}), params = extra_args{k}; end
    end
  end

  J_prev = []; Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
  h = min(max_step, 1e-2);
  LU_cached = struct(); LU_cached.ilu_fallback = false;

  while t < tspan(2) - 1e-12
    h = min(h, tspan(2) - t);
    K = zeros(n_state,3); Xstage = repmat(x,1,3);
    for stage = 1:3
      xi = x;
      for newt = 1:8
        ti = t + sum(A(stage,1:stage))*h;
        x_tau = HistoryCache_eval(H, ti - h_delay);
        fval = f_rhs(ti, xi, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
        R = xi - x - h*(A(stage,stage)*fval + sum(A(stage,1:stage-1).*K(:,1:stage-1),2));
        J = analytic_jacobian_reactor(xi, params, Lmat);
        M = speye(n_state) - h*A(stage,stage)*J;

        % LU caching and robust solve
        recomputeLU = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
        if recomputeLU
          try
            [Lfac,Ufac,Pfac,Qfac,Rfac] = lu(M);
            LU_cached.L = Lfac; LU_cached.U = Ufac; LU_cached.P = Pfac;
            LU_cached.Q = Qfac; LU_cached.R = Rfac; LU_cached.ilu_fallback = false;
            J_prev = J;
          catch
            warning('Sparse LU failed; using ILU+GMRES fallback.');
            setup.type = 'ilutp'; setup.droptol = 1e-3;
            [Lilu,Uilu] = ilu(M, setup);
            LU_cached.L = Lilu; LU_cached.U = Uilu; LU_cached.P = []; LU_cached.Q = []; LU_cached.R = [];
            LU_cached.ilu_fallback = true; J_prev = J;
          end
        end

        RHS = -R;
        if LU_cached.ilu_fallback
          tol_lin = 1e-6; maxit_lin = 200;
          try
            [delta, flag] = gmres(M, RHS, [], tol_lin, maxit_lin, LU_cached.L, LU_cached.U);
            if flag ~= 0
              delta = M \ RHS;
            end
          catch
            delta = M \ RHS;
          end
        else
          if ~isempty(LU_cached.R)
            y = LU_cached.P * RHS;
            z = LU_cached.L \ y;
            w = LU_cached.U \ z;
            w = LU_cached.R \ w;
            delta = LU_cached.Q * w;
          else
            y = LU_cached.P * RHS;
            z = LU_cached.L \ y;
            w = LU_cached.U \ z;
            delta = LU_cached.Q * w;
          end
        end

        xi = xi + delta; stats.newton_iters = stats.newton_iters + 1; stats.linear_solves = stats.linear_solves + 1;
        if norm(delta, inf) < 1e-6, break; end
      end
      K(:,stage) = fval; Xstage(:,stage) = xi;
    end
    x_new = x + h*(K * b);
    err_est = norm(x_new - Xstage(:,3), inf);
    tol = atol + rtol * max(norm(x,inf), norm(x_new,inf));
    if err_est <= tol
      t = t + h; x = x_new; T_out(end+1,1) = t; X_out(end+1,:) = x.'; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
      if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
      h = min(max_step, max(1e-8, fac*h));
    else
      h = max(1e-8, 0.5*h);
    end
  end
end

function [T_out, X_out, stats] = simulate_dde_rosenbrock(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts)
  if nargin < 7, opts = struct(); end
  atol = getfielddef(opts,'atol',1e-5);
  rtol = getfielddef(opts,'rtol',1e-5);
  max_step = getfielddef(opts,'max_step',(tspan(2)-tspan(1))/100);
  interp_grid = getfielddef(opts,'interp_grid',0.02);

  gamma = 0.5;
  h_eff = max(h_delay, 1e-8);
  hist_times = (-h_eff:interp_grid:0);
  hist_states = zeros(numel(hist_times), n_state);
  for k = 1:numel(hist_times), hist_states(k,:) = x0_fun(hist_times(k)).'; end
  H = HistoryCache_create(hist_times, hist_states, 200000);

  t = tspan(1); x = x0_fun(0); x = x(:);
  T_out = t; X_out = x.'; stats.rhs_evals = 0; stats.linear_solves = 0; stats.steps = 0;

  Lmat = []; params = [];
  if ~isempty(extra_args)
    for k = 1:numel(extra_args)
      if issparse(extra_args{k}) && isempty(Lmat), Lmat = extra_args{k};
      elseif isstruct(extra_args{k}), params = extra_args{k}; end
    end
  end

  J_prev = []; Jtol = getfielddef(opts,'jacobian_reuse_tol',1e-8);
  h = min(max_step, 1e-2);
  LU_cached = struct(); LU_cached.ilu_fallback = false;

  while t < tspan(2) - 1e-12
    h = min(h, tspan(2) - t);
    x_tau = HistoryCache_eval(H, t - h_delay);
    J = analytic_jacobian_reactor(x, params, Lmat);
    M = speye(n_state) - gamma*h*J;

    recomputeLU = isempty(J_prev) || norm(J - J_prev, inf) > Jtol;
    if recomputeLU
      try
        [Lfac,Ufac,Pfac,Qfac,Rfac] = lu(M);
        LU_cached.L = Lfac; LU_cached.U = Ufac; LU_cached.P = Pfac;
        LU_cached.Q = Qfac; LU_cached.R = Rfac; LU_cached.ilu_fallback = false;
        J_prev = J;
      catch
        warning('Sparse LU failed in Rosenbrock; using ILU+GMRES fallback.');
        setup.type = 'ilutp'; setup.droptol = 1e-3;
        [Lilu,Uilu] = ilu(M, setup);
        LU_cached.L = Lilu; LU_cached.U = Uilu; LU_cached.P = []; LU_cached.Q = []; LU_cached.R = [];
        LU_cached.ilu_fallback = true; J_prev = J;
      end
    end

    f1 = f_rhs(t, x, x_tau); stats.rhs_evals = stats.rhs_evals + 1;
    if LU_cached.ilu_fallback
      k1 = gmres(M, f1, [], 1e-8, 200, LU_cached.L, LU_cached.U);
    else
      y = LU_cached.P * f1;
      z = LU_cached.L \ y;
      w = LU_cached.U \ z;
      if ~isempty(LU_cached.R), w = LU_cached.R \ w; end
      k1 = LU_cached.Q * w;
    end
    t2 = t + 0.5*h; x2 = x + 0.5*h*k1;
    x2_tau = HistoryCache_eval(H, t2 - h_delay);
    f2 = f_rhs(t2, x2, x2_tau); stats.rhs_evals = stats.rhs_evals + 1;
    rhs2 = f2 + (1/(2*h))*k1;
    if LU_cached.ilu_fallback
      k2 = gmres(M, rhs2, [], 1e-8, 200, LU_cached.L, LU_cached.U);
    else
      y2 = LU_cached.P * rhs2;
      z2 = LU_cached.L \ y2;
      w2 = LU_cached.U \ z2;
      if ~isempty(LU_cached.R), w2 = LU_cached.R \ w2; end
      k2 = LU_cached.Q * w2;
    end
    x_new = x + h*(k1 + k2)/2;
    err_est = norm(x_new - (x + h*k1), inf);
    tol = atol + rtol * max(norm(x,inf), norm(x_new,inf));
    if err_est <= tol
      t = t + h; x = x_new; T_out(end+1,1) = t; X_out(end+1,:) = x.'; H = HistoryCache_append(H, t, x); stats.steps = stats.steps + 1;
      if err_est == 0, fac = 2.0; else fac = 0.9*(tol/err_est)^(1/3); end
      h = min(max_step, max(1e-8, fac*h));
    else
      h = max(1e-8, 0.5*h);
    end
  end
end

%% ------------------------------------------------------------------------
%% 3D plotting utilities: optimized to avoid heavy isosurface on large grids
%% ------------------------------------------------------------------------
function create_3d_plots(P, Lmat, Tvec, Xmat, title_prefix)
  Nz = P.grid.Nz; Nr = P.grid.Nr; N = P.Ncells;
  z = linspace(0, getfielddef(P,'geometry.L',getfielddef(P,'geometry.tube_length',1.0)), Nz);
  r = linspace(0, getfielddef(P,'geometry.R',getfielddef(P,'geometry.tube_radius',1.0)), Nr);
  c_final = Xmat(end, 1:N)'; T_final = Xmat(end, N+1:end)';
  C = reshape(c_final, [Nr, Nz]); Tm = reshape(T_final, [Nr, Nz]);

  max_cells_for_isosurface = 200000;
  if P.Ncells <= max_cells_for_isosurface
    nang = 36;
    theta = linspace(0,2*pi,nang);
    [TH,Zg] = meshgrid(theta, z);
    Xcart = zeros(Nr, nang, Nz); Ycart = Xcart; Zcart = zeros(Nr, nang, Nz);
    Ccart = zeros(Nr, nang, Nz); Tcart = zeros(Nr, nang, Nz);
    for iz = 1:Nz
      for ir = 1:Nr
        for it = 1:nang
          Xcart(ir,it,iz) = r(ir)*cos(theta(it));
          Ycart(ir,it,iz) = r(ir)*sin(theta(it));
          Zcart(ir,it,iz) = z(iz);
          Ccart(ir,it,iz) = C(ir,iz);
          Tcart(ir,it,iz) = Tm(ir,iz);
        end
      end
    end
    figure('Name',[title_prefix ' 3D Temperature isosurface']);
    p = patch(isosurface(Xcart, Ycart, Zcart, Tcart, mean(Tcart(:)))); isonormals(Xcart,Ycart,Zcart,Tcart,p);
    set(p,'FaceColor','red','EdgeColor','none','FaceAlpha',0.7);
    camlight; lighting gouraud; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
    title([title_prefix ' — Temperature isosurface (mean)']);

    figure('Name',[title_prefix ' 3D Concentration isosurface']);
    p2 = patch(isosurface(Xcart, Ycart, Zcart, Ccart, mean(Ccart(:))));
    set(p2,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.6);
    camlight; lighting gouraud; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
    title([title_prefix ' — Concentration isosurface (mean)']);
  else
    figure('Name',[title_prefix ' Concentration axial slice']);
    imagesc(z, r, C); colorbar; axis xy; xlabel('z (m)'); ylabel('r (m)');
    title([title_prefix ' — Concentration axial slice']);
    figure('Name',[title_prefix ' Centerline Temperature']);
    center_idx = ceil(Nr/2);
    plot(z, Tm(center_idx,:)); xlabel('z (m)'); ylabel('T (K)');
    title([title_prefix ' — Centerline Temperature']);
  end

  pressure_proxy = compute_pressure_proxy(P, C, z, r);
  figure('Name',[title_prefix ' Pressure proxy (axial)']);
  imagesc(z, r, pressure_proxy); colorbar; axis xy; xlabel('z (m)'); ylabel('r (m)');
  title([title_prefix ' — Pressure proxy (axial)']);
end

function p = compute_pressure_proxy(P, C, z, r)
  Nz = numel(z); Nr = numel(r);
  Kbase = 1e3;
  p = zeros(Nr, Nz);
  for iz = 1:Nz
    for ir = 1:Nr
      p(ir,iz) = Kbase * (1 - getfielddef(P,'phys.porosity',0.4)) * (1 + C(ir,iz));
    end
  end
end

%% ------------------------------------------------------------------------
%% Utility functions
%% ------------------------------------------------------------------------
function v = getfielddef(s,f,d)
  if isfield(s,f), v = s.(f); else v = d; end
end

function err = rel_error_vs_ref(Tm, Xm, Tref, Xref)
  Xinterp = interp1(Tm, Xm, Tref);
  err = norm(Xinterp(:) - Xref(:)) / max(1e-12, norm(Xref(:)));
end



function industrial_reactors_highres_full_fixed()

  %% ---------------- User setters (geometry, industrial params, simulation time) ----------------
  A = set_packed_bed_geometry(12.0, 1.6, 240, 60);   % L (m), R (m), Nz, Nr
  B = set_multi_tubular_geometry(18.0, 0.15, 400, 32); % tube_length (m), tube_radius (m), Nz, Nr

  A = set_packed_bed_industrial_params(A, ...
       struct('Dc',1.0e-5,'DT',5.0e-6,'k1',10.0,'Ea',6.0e4,'k2',0.5,'beta',5.0,'gamma',1.0,'porosity',0.4));
  B = set_multi_tubular_industrial_params(B, ...
       struct('Dc',5.0e-6,'DT',2.0e-6,'k1',1.2,'Ea',4.5e4,'k2',0.15,'beta',2.0,'gamma',0.6,'porosity',0.45));

  A.Tfinal = 60.0;
  B.Tfinal = 120.0;

  A.solver = default_solver_settings('rosenbrock',1e-7,1e-7,0.02,0.005,1e-8);
  B.solver = default_solver_settings('sdirk',1e-7,1e-7,0.02,0.005,1e-8);

  %% ---------------- Build spatial operators and initial conditions ----------------
  L_A = build_axisymmetric_laplacian(A.grid.Nz, A.grid.Nr, A.dx, A.dr);
  L_B = build_axisymmetric_laplacian(B.grid.Nz, B.grid.Nr, B.dx, B.dr);

  x0A_fun = @(t) [A.ic.c0; A.ic.T0];
  x0B_fun = @(t) [B.ic.c0; B.ic.T0];

  %% ---------------- Run high-resolution simulations and benchmarks ----------------
  fprintf('Running high-resolution simulations (this may take time)...\n');

  tspanA = [0, A.Tfinal];
  opts_refA = struct('atol',1e-9,'rtol',1e-9,'max_step',A.Tfinal/1000,'interp_grid',A.history.grid,'jacobian_reuse_tol',A.solver.jacobian_reuse_tol);
  [TrefA, XrefA, ~] = simulate_dde_sdirk(@(t,x,xtau) rhs_packed(t,x,xtau,L_A,A), x0A_fun, A.delay.h, tspanA, 2*A.Ncells, {L_A,A}, opts_refA);

  opts_rosA = A.solver.opts; opts_rosA.jacobian_reuse_tol = A.solver.jacobian_reuse_tol;
  tic; [TrosA, XrosA, statsRosA] = simulate_dde_rosenbrock(@(t,x,xtau) rhs_packed(t,x,xtau,L_A,A), x0A_fun, A.delay.h, tspanA, 2*A.Ncells, {L_A,A}, opts_rosA); timeRosA = toc;

  opts_sdirkA = A.solver.opts; opts_sdirkA.jacobian_reuse_tol = A.solver.jacobian_reuse_tol;
  tic; [TsdirkA, XsdirkA, statsSdirkA] = simulate_dde_sdirk(@(t,x,xtau) rhs_packed(t,x,xtau,L_A,A), x0A_fun, A.delay.h, tspanA, 2*A.Ncells, {L_A,A}, opts_sdirkA); timeSdirkA = toc;

  opts_rk4A = struct('atol',1e-6,'rtol',1e-6,'max_step',A.Tfinal/200,'interp_grid',A.history.grid);
  tic; [Trk4A, Xrk4A, statsRk4A] = simulate_dde_rk4_adaptive(@(t,x,xtau) rhs_packed(t,x,xtau,L_A,A), x0A_fun, A.delay.h, tspanA, 2*A.Ncells, {L_A,A}, opts_rk4A); timeRk4A = toc;

  errRosA   = rel_error_vs_ref(TrosA, XrosA, TrefA, XrefA);
  errSdirkA = rel_error_vs_ref(TsdirkA, XsdirkA, TrefA, XrefA);
  errRk4A   = rel_error_vs_ref(Trk4A, Xrk4A, TrefA, XrefA);

  tspanB = [0, B.Tfinal];
  opts_refB = struct('atol',1e-9,'rtol',1e-9,'max_step',B.Tfinal/1000,'interp_grid',B.history.grid,'jacobian_reuse_tol',B.solver.jacobian_reuse_tol);
  [TrefB, XrefB, ~] = simulate_dde_sdirk(@(t,x,xtau) rhs_packed(t,x,xtau,L_B,B), x0B_fun, B.delay.h, tspanB, 2*B.Ncells, {L_B,B}, opts_refB);

  opts_sdirkB = B.solver.opts; opts_sdirkB.jacobian_reuse_tol = B.solver.jacobian_reuse_tol;
  tic; [TsdirkB, XsdirkB, statsSdirkB] = simulate_dde_sdirk(@(t,x,xtau) rhs_packed(t,x,xtau,L_B,B), x0B_fun, B.delay.h, tspanB, 2*B.Ncells, {L_B,B}, opts_sdirkB); timeSdirkB = toc;

  opts_rosB = B.solver.opts; opts_rosB.jacobian_reuse_tol = B.solver.jacobian_reuse_tol;
  tic; [TrosB, XrosB, statsRosB] = simulate_dde_rosenbrock(@(t,x,xtau) rhs_packed(t,x,xtau,L_B,B), x0B_fun, B.delay.h, tspanB, 2*B.Ncells, {L_B,B}, opts_rosB); timeRosB = toc;

  opts_rk4B = struct('atol',1e-6,'rtol',1e-6,'max_step',B.Tfinal/200,'interp_grid',B.history.grid);
  tic; [Trk4B, Xrk4B, statsRk4B] = simulate_dde_rk4_adaptive(@(t,x,xtau) rhs_packed(t,x,xtau,L_B,B), x0B_fun, B.delay.h, tspanB, 2*B.Ncells, {L_B,B}, opts_rk4B); timeRk4B = toc;

  errSdirkB = rel_error_vs_ref(TsdirkB, XsdirkB, TrefB, XrefB);
  errRosB   = rel_error_vs_ref(TrosB, XrosB, TrefB, XrefB);
  errRk4B   = rel_error_vs_ref(Trk4B, Xrk4B, TrefB, XrefB);

  %% ---------------- 3D visualization: reconstruct Cartesian mesh and plot ----------------
  create_3d_plots(A, L_A, TrosA, XrosA, 'Packed‑bed gaseous fuel reactor (A)');
  create_3d_plots(B, L_B, TsdirkB, XsdirkB, 'Multi‑tubular polymer reactor (B)');

  %% ---------------- Benchmark summary (console + save) ----------------
  fprintf('\n=== Benchmark summary (high‑res) ===\n');
  fprintf('Reactor A (Packed‑bed): RK4 time=%.2fs err=%.3e | SDIRK time=%.2fs err=%.3e | Rosenbrock time=%.2fs err=%.3e\n', ...
          timeRk4A, errRk4A, timeSdirkA, errSdirkA, timeRosA, errRosA);
  fprintf('Reactor B (Multi‑tubular): RK4 time=%.2fs err=%.3e | SDIRK time=%.2fs err=%.3e | Rosenbrock time=%.2fs err=%.3e\n', ...
          timeRk4B, errRk4B, timeSdirkB, errSdirkB, timeRosB, errRosB);

  save('highres_results_full_fixed.mat','A','B','TrefA','XrefA','TrefB','XrefB','-v7.3');

  fprintf('All done. Results saved to highres_results_full_fixed.mat\n');
end

industrial_reactors_highres_full_fixed()


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
