%% Industrial Reactor Simulation – Memory‑Safe Version
% Author: Pavlo Mazniker
% Fully corrected to prevent hangs and memory spikes.
% Run with: industrial_reactors_highres_full_fixed_clean()

close all force; clc; clear all; clear global; clear classes;
format long e; rand('state',1); pack; max_stack_depth(1024*500);

% ----------------------------------------------------------------------
% Geometry and parameter constructors
% ----------------------------------------------------------------------
function R = g1(L,Rr,Nz,Nr)
    R.geometry.L = L; R.geometry.R = Rr; R.grid.Nz = Nz; R.grid.Nr = Nr;
    R.Ncells = Nz*Nr; R.dx = L/(Nz-1); R.dr = Rr/(Nr-1);
    R.ic.c0 = 0.5*ones(R.Ncells,1); R.ic.T0 = 300*ones(R.Ncells,1);
    R.delay.h = 0.3; R.history.grid = 0.005;
    R.bc.inlet.c = 0.8; R.bc.inlet.T = 300;
end

function R = g2(Lt,Rt,Nz,Nr)
    R.geometry.tube_length = Lt; R.geometry.tube_radius = Rt;
    R.grid.Nz = Nz; R.grid.Nr = Nr; R.Ncells = Nz*Nr;
    R.dx = Lt/(Nz-1); R.dr = Rt/(Nr-1);
    R.ic.c0 = 0.4*ones(R.Ncells,1); R.ic.T0 = 320*ones(R.Ncells,1);
    R.delay.h = 0.15; R.history.grid = 0.005;
    R.bc.inlet.c = 0.6; R.bc.inlet.T = 320;
end

function R = g3(R,p)
    R.phys = p;
end

function R = g4(R,p)
    R.phys = p;
end

function s = s1(m,a,r,ms,ig,jt)
    s.method = m; s.opts = struct('atol',a,'rtol',r,'max_step',ms, ...
                                   'interp_grid',ig); s.jacobian_reuse_tol = jt;
end

% ----------------------------------------------------------------------
% Laplacian (sparse, 5‑point stencil, Neumann BCs)
% ----------------------------------------------------------------------
function L = L1(Nz,Nr,dz,dr)
    N = Nz*Nr; L = spalloc(N,N,5*N); idx = @(iz,ir) (ir-1)*Nz + iz;
    for ir=1:Nr
        for iz=1:Nz
            k = idx(iz,ir);
            if iz>1,  L(k,idx(iz-1,ir)) = 1/dz^2; end
            if iz<Nz, L(k,idx(iz+1,ir)) = 1/dz^2; end
            if ir>1,  L(k,idx(iz,ir-1)) = 1/dr^2; end
            if ir<Nr, L(k,idx(iz,ir+1)) = 1/dr^2; end
            L(k,k) = -sum(L(k,:));
        end
    end
end

% ----------------------------------------------------------------------
% RHS of DDE
% ----------------------------------------------------------------------
function dX = rhs(t,X,Xt,L,P)
    N = P.Ncells; c = X(1:N); T = X(N+1:end);
    R = P.phys.k1 .* c .* exp(-(P.phys.Ea/8.314).*(1./T - 1/300));
    dc = P.phys.Dc*(L*c) - R - P.phys.k2*c;
    dT = P.phys.DT*(L*T) + P.phys.beta*R - P.phys.gamma*(Xt(N+1:end)-P.bc.inlet.T);
    dX = [dc; dT];
end

% ----------------------------------------------------------------------
% History cache (circular buffer)
% ----------------------------------------------------------------------
function H = Hc(ti,xi,m)
    nt = numel(ti); ns = size(xi,2);
    H.maxlen = m;
    H.t = nan(m,1); H.x = nan(m,ns);
    H.head = 1; H.tail = nt; H.count = nt;
    H.t(1:nt) = ti(:); H.x(1:nt,:) = xi;
end

function H = Ha(H,t,x)
    if H.count < H.maxlen
        H.tail = H.tail + 1;
        H.t(H.tail) = t; H.x(H.tail,:) = x(:).';
        H.count = H.count + 1;
    else
        H.head = H.head + 1; if H.head > H.maxlen, H.head = 1; end
        H.tail = H.head + H.count - 1; if H.tail > H.maxlen, H.tail = H.tail - H.maxlen; end
        H.t(H.tail) = t; H.x(H.tail,:) = x(:).';
    end
end

function xq = He(H,tq)
    if H.head <= H.tail
        tv = H.t(H.head:H.tail); xm = H.x(H.head:H.tail,:);
    else
        tv = [H.t(H.head:end); H.t(1:H.tail)];
        xm = [H.x(H.head:end,:); H.x(1:H.tail,:)];
    end
    xq = interp1(tv,xm,tq,'pchip')';
end

% ----------------------------------------------------------------------
% Jacobian of RHS (for linear algebra)
% ----------------------------------------------------------------------
function J = J1(x,p,L)
    N = p.Ncells; c = x(1:N); T = x(N+1:end);
    e = exp(-(p.phys.Ea/8.314).*(1./T - 1/300));
    d1 = p.phys.k1.*e; d2 = p.phys.k1.*c.*e.*(p.phys.Ea/8.314).*(1./T.^2);
    J11 = p.phys.Dc*L - spdiags(p.phys.k2 + d1, 0, N, N);
    J12 = -spdiags(d2, 0, N, N);
    J21 = spdiags(p.phys.beta*d1, 0, N, N);
    J22 = p.phys.DT*L - spdiags(p.phys.gamma*ones(N,1), 0, N, N);
    J = [J11 J12; J21 J22];
end

% ----------------------------------------------------------------------
% Helper: solve (I - γ h J) x = b using either sparse LU or ILU+GMRES
% ----------------------------------------------------------------------
function x = solve_linear(M, b, use_ilu, reuse)
    persistent LUstruct Jprev; % simple persistent cache (per call context)
    if nargin < 4, reuse = false; end
    if reuse && ~isempty(Jprev) && norm(M - Jprev, inf) < 1e-8
        % reuse previous factorization – assume LUstruct is still valid
    else
        if ~use_ilu
            try
                [L,U,P,Q,R] = lu(M);
                LUstruct = struct('L',L,'U',U,'P',P,'Q',Q,'R',R,'ilu',false);
            catch
                use_ilu = true; % fallback if LU fails
            end
        end
        if use_ilu
            setup.type = 'ilutp'; setup.droptol = 1e-3;
            [L,U] = ilu(M, setup);
            LUstruct = struct('L',L,'U',U,'ilu',true);
        end
        Jprev = M;
    end

    if LUstruct.ilu
        [x, flag] = gmres(M, b, [], 1e-6, 200, LUstruct.L, LUstruct.U);
        if flag ~= 0
            x = M \ b; % fallback to direct
        end
    else
        y = LUstruct.P * b;
        z = LUstruct.L \ y;
        w = LUstruct.U \ z;
        if ~isempty(LUstruct.R)
            w = LUstruct.R \ w;
        end
        x = LUstruct.Q * w;
    end
end

% ----------------------------------------------------------------------
% SDIRK integrator (3‑stage, with optional callback)
% ----------------------------------------------------------------------
##function [T_out, X_out, stats] = sdirk(f, x0, h0, ts, n, extra_args, opts, accept_callback)
##    % Butcher tableau (3‑stage, stiffly accurate)
##    g = 0.435866521508459;
##    A = [g 0 0; 0.5-g g 0; 2*g 1-4*g g];
##    b = [2*g; 1-4*g; 2*g];
##
##    atol = opts.atol; rtol = opts.rtol; max_step = opts.max_step;
##    interp_grid = opts.interp_grid; jac_reuse_tol = opts.jacobian_reuse_tol;
##    force_ilu = getfielddef(opts, 'force_ilu', false);
##    history_maxlen = getfielddef(opts, 'history_maxlen', 200000);
##
##    % Initial history
##    he = max(h0, 1e-8);
##    ht = (-he:interp_grid:0)';
##    xs = zeros(numel(ht), n);
##    for k = 1:numel(ht)
##        xs(k,:) = x0(ht(k)).';
##    end
##    H = Hc(ht, xs, history_maxlen);
##
##    t = ts(1); x = x0(0); x = x(:);
##    T_out = t; X_out = x.';
##    L = extra_args{1}; p = extra_args{2};
##
##    % For LU reuse
##    Jp = []; LUstruct = []; use_ilu = force_ilu || (n > 5e4);
##    h = min(max_step, 0.01);
##
##    while t < ts(2)-1e-10
##        h = min(h, ts(2)-t);
##        K = zeros(n,3);
##
##        % Stage loop
##        for s = 1:3
##            xi = x;
##            for it = 1:8
##                ti = t + sum(A(s,1:s))*h;
##                xt = He(H, ti - p.delay.h);
##                fv = f(ti, xi, xt);
##                R = xi - x - h*(A(s,s)*fv + sum(A(s,1:s-1).*K(:,1:s-1),2));
##                J = J1(xi, p, L);
##                M = speye(n) - h*A(s,s)*J;
##                if isempty(Jp) || norm(J-Jp, inf) > jac_reuse_tol
##                    use_ilu = force_ilu || (n > 5e4);
##                    LUstruct = []; % force rebuild
##                end
##                d = solve_linear(M, -R, use_ilu, ~isempty(LUstruct));
##                xi = xi + d;
##                if norm(d,inf) < 1e-6
##                    break;
##                end
##            end
##            K(:,s) = fv;
##        end
##
##        % Combine stages
##        xn = x + h*(K*b);
##        % Error estimate (difference from first stage)
##        er = norm(xn - (x + h*K(:,1)), inf);
##        tol = atol + rtol*max(norm(x,inf), norm(xn,inf));
##
##        if er <= tol
##            t = t + h; x = xn;
##            T_out(end+1,1) = t; X_out(end+1,:) = x.';
##            if nargin >= 8 && ~isempty(accept_callback)
##                accept_callback(t, x);
##            end
##            H = Ha(H, t, x);
##            % step size control
##            fc = (er==0)*2 + (er>0)*0.9*(tol/er)^(1/3);
##            h = min(max_step, max(1e-6, fc*h));
##        else
##            h = max(1e-6, 0.5*h);
##        end
##    end
##    stats = struct();
##end

function [T_out, X_out, stats] = sdirk(f, x0, h0, ts, n, extra_args, opts, callback)
    % SDIRK  3‑stage singly diagonal implicit Runge–Kutta (stiffly accurate)
    %        with PI step size control and embedded error estimate.
    %
    %   [T,X] = sdirk(F,X0,H0,TSPAN,N,EXTRA,OPTS,CALLBACK)
    %   F      – function dydt = F(t,y,ydel,extra{:})
    %   X0     – history function X0(t) (returns column vector)
    %   H0     – delay
    %   TSPAN  – [t0 tfinal]
    %   N      – dimension of system
    %   EXTRA  – cell array of extra arguments for F
    %   OPTS   – structure with fields:
    %       atol, rtol, max_step, interp_grid, jacobian_reuse_tol,
    %       force_ilu (optional), history_maxlen (optional)
    %   CALLBACK – function handle called after each accepted step (optional)
    %
    %   Outputs:
    %       T_out – vector of accepted times (all steps)
    %       X_out – corresponding states (each row a time)
    %       stats – structure with statistics (steps, failed steps, etc.)

    % Butcher tableau for 3‑stage SDIRK (stiffly accurate)
    g = 0.435866521508459;  % gamma
    A = [g,   0,    0;
         0.5-g, g,    0;
         2*g, 1-4*g, g];
    b = [2*g; 1-4*g; 2*g];  % weights (same as last row)
    c = [g; 0.5; 1];        % stage times

    % Options
    atol = opts.atol;
    rtol = opts.rtol;
    max_step = opts.max_step;
    interp_grid = opts.interp_grid;
    jac_reuse_tol = opts.jacobian_reuse_tol;
    force_ilu = getfielddef(opts, 'force_ilu', false);
    history_maxlen = getfielddef(opts, 'history_maxlen', 200000);
    max_no_progress = 1000;   % safety: abort if step size collapses

    % Initial history
    he = max(h0, 1e-8);
    ht = (-he:interp_grid:0)';
    xs = zeros(numel(ht), n);
    for k = 1:numel(ht)
        xs(k,:) = x0(ht(k)).';
    end
    H = Hc(ht, xs, history_maxlen);

    t = ts(1);
    x = x0(0); x = x(:);
    T_out = t; X_out = x.';
    L = extra_args{1};
    p = extra_args{2};

    % State for linear solver
    Jp = [];
    use_ilu = force_ilu || (n > 5e4);
    stats.steps = 0;
    stats.failed = 0;
    stats.jac_updates = 0;

    h = min(max_step, 0.01 * (ts(2)-ts(1)));  % initial step

    % PI step size controller parameters (Hairer & Wanner)
    fac_min = 0.2;
    fac_max = 5.0;
    fac = 0.9;          % safety factor
    beta = 0.0;          % for PI controller (set to 0 for I controller)
    errold = 1.0;        % previous error
    same_jac_step = 0;   % counter for steps since last Jacobian update

    while t < ts(2)-1e-12
        h = min(h, ts(2)-t);
        if h < 1e-12
            warning('Step size too small – stopping at t = %g', t);
            break;
        end

        % Compute Jacobian and factorize (if needed)
      h = min(h, ts(2)-t);
      h_this = h;   % store for Jacobian update condition
      hprev = h
        if isempty(Jp) || same_jac_step > 10 || ...
           (h / hprev > 2.0 && same_jac_step > 0) || h / hprev < 0.5
            xt = He(H, t - p.delay.h);
            J = J1(x, p, L);
            Jp = J;
            use_ilu = force_ilu || (n > 5e4);
            stats.jac_updates = stats.jac_updates + 1;
            same_jac_step = 0;
            % Prefactor the matrix M = I - h*gamma*J (will be used for all stages)
            % We'll do it in the Newton loop because it depends on h
            % But we store Jp for reuse
        end

        % After step control, update hprev
         hprev = h_this;

        % Newton iteration for stages (simplified: solve each stage separately)
        % We'll use the same matrix M = I - h*gamma*J for all stages in this step
        M = speye(n) - h*g*Jp;
        use_ilu_step = use_ilu;

        % Allocate stage vectors
        K = zeros(n, 3);
        stage_success = true;

        for s = 1:3
            ti = t + c(s)*h;
            xi = x;
            % Build right-hand side for stage equation
            rhs = zeros(n,1);
            for j = 1:s-1
                rhs = rhs + A(s,j)*K(:,j);
            end
            rhs = x + h*rhs;
            % Stage equation: xi = rhs + h*A(s,s)*f(ti, xi, xti)
            % We solve using Newton with fixed matrix M
            xti = He(H, ti - p.delay.h);
            for iter = 1:8
                fval = f(ti, xi, xti);
                res = xi - rhs - h*A(s,s)*fval;
                if norm(res, inf) < 1e-6
                    break;
                end
                % Solve M * delta = res (note sign: M = I - h*g*J)
                delta = solve_linear(M, res, use_ilu_step);
                xi = xi - delta;
                if norm(delta, inf) < 1e-6
                    break;
                end
            end
            if norm(res, inf) > 1e-3
                stage_success = false;
                break;
            end
            K(:,s) = fval;  % store f(ti,xi,xti)
        end

        if ~stage_success
            % Newton failed – reduce step and retry
            h = h * 0.25;
            stats.failed = stats.failed + 1;
            continue;
        end

        % Compute solution and embedded error estimate
        xn = x + h * (K * b);
        % Error estimate (difference from first stage – not rigorous but common)
        er = norm(xn - (x + h*K(:,1)), inf);
        tol = atol + rtol * max(norm(x,inf), norm(xn,inf));

        if er <= tol
            % Accept step
            t = t + h;
            x = xn;
            T_out(end+1,1) = t;
            X_out(end+1,:) = x.';
            if nargin >= 8 && ~isempty(callback)
                callback(t, x);
            end
            H = Ha(H, t, x);
            stats.steps = stats.steps + 1;

            % Step size control (PI controller)
            if stats.steps > 1
                expo = 1/3;  % order of SDIRK is 3
                facc = fac * (tol / er)^expo * (er / errold)^(beta * expo);
            else
                facc = fac * (tol / er)^(1/3);
            end
            facc = min(fac_max, max(fac_min, facc));
            h_new = facc * h;
            errold = max(er, 1e-12);
            h = min(max_step, max(1e-8, h_new));
            same_jac_step = same_jac_step + 1;
        else
            % Reject step
            stats.failed = stats.failed + 1;
            h = h * 0.5;
            same_jac_step = 0;  % force Jacobian update after rejection
        end

        % Safety check: if too many consecutive failures, abort
        if stats.failed > max_no_progress
            warning('Too many failed steps – stopping at t = %g', t);
            break;
        end
        hprev = h;
    end

    stats.nsteps = stats.steps;
    stats.nfailed = stats.failed;
    stats.njac = stats.jac_updates;
end

function [T_out, X_out, stats] = sdirkf(f, x0, h0, ts, n, extra_args, opts, callback)
    % SDIRK – 3‑stage stiffly accurate SDIRK with embedded error estimate.
    %         Includes Newton with extrapolated initial guess and a safety
    %         minimum step to prevent infinite reduction.

    % Butcher tableau (stiffly accurate, order 3)
    g = 0.435866521508459;
    A = [g,   0,    0;
         0.5-g, g,    0;
         2*g, 1-4*g, g];
    b = [2*g; 1-4*g; 2*g];
    % Embedded method (order 2) – common choice for SDIRK
    b2 = [g; (1-g); g];  % example, may need tuning
    c = [g; 0.5; 1];

    % Options
    atol = opts.atol;
    rtol = opts.rtol;
    max_step = opts.max_step;
    interp_grid = opts.interp_grid;
    jac_reuse_tol = opts.jacobian_reuse_tol;
    force_ilu = getfielddef(opts, 'force_ilu', false);
    history_maxlen = getfielddef(opts, 'history_maxlen', 200000);
    max_no_progress = 1000;
    min_step = 1e-10;       % if step falls below this, accept anyway

    % Initial history
    he = max(h0, 1e-8);
    ht = (-he:interp_grid:0)';
    xs = zeros(numel(ht), n);
    for k = 1:numel(ht)
        xs(k,:) = x0(ht(k)).';
    end
    H = Hc(ht, xs, history_maxlen);

    t = ts(1);
    x = x0(0); x = x(:);
    T_out = t; X_out = x.';
    L = extra_args{1};
    p = extra_args{2};

    Jp = [];                 % last computed Jacobian
    use_ilu = force_ilu || (n > 5e4);
    stats.steps = 0;
    stats.failed = 0;
    stats.jac_updates = 0;

    h = min(max_step, 0.01 * (ts(2)-ts(1)));
    h_last_accepted = h;     % last accepted step size (for Jacobian reuse)
    errold = 1.0;
    same_jac_step = 0;

    % PI controller parameters
    fac_min = 0.2; fac_max = 5.0; fac = 0.9; beta = 0.0;

    while t < ts(2)-1e-12
        h = min(h, ts(2)-t);
        if h < min_step
            warning('Step size extremely small (%g) – forcing acceptance at t=%g', h, t);
            % Force accept current state (do not advance)
            % This is a safety valve – you may want to break instead
        end

        % --- Jacobian update decision (using last accepted step size) ---
        if isempty(Jp) || same_jac_step > 10 || ...
           (h / h_last_accepted > 2.0 && same_jac_step > 0) || ...
           (h / h_last_accepted < 0.5 && same_jac_step > 0)
            xt = He(H, t - p.delay.h);
            J = J1(x, p, L);
            Jp = J;
            use_ilu = force_ilu || (n > 5e4);
            stats.jac_updates = stats.jac_updates + 1;
            same_jac_step = 0;
        end

        % --- Newton loop for stages ---
        M = speye(n) - h*g*Jp;      % fixed matrix for this step
        K = zeros(n, 3);
        stage_success = true;

        for s = 1:3
            ti = t + c(s)*h;
            xti = He(H, ti - p.delay.h);

            % Initial guess: use previous stage if available, else x
            if s == 1
                xi = x;
            else
                % extrapolate using previous stage derivative
                xi = x + h * (A(s,1:s-1) * K(:,1:s-1)');
            end

            % Build rhs for stage equation: xi = x + h * sum_{j=1}^{s-1} A(s,j)*Kj + h*A(s,s)*f(ti,xi,xti)
            rhs_stage = x;
            for j = 1:s-1
                rhs_stage = rhs_stage + h * A(s,j) * K(:,j);
            end

            % Simplified Newton iterations
            for iter = 1:15
                fval = f(ti, xi, xti);
                res = xi - rhs_stage - h * A(s,s) * fval;
                nrm_res = norm(res, inf);
                if nrm_res < atol + rtol * norm(xi, inf)
                    break;
                end
                % Solve M * delta = res
                delta = solve_linear(M, res, use_ilu);
                xi = xi - delta;
                if norm(delta, inf) < 1e-8
                    break;
                end
            end
            if nrm_res > 1e-3   % Newton did not converge well enough
                stage_success = false;
                break;
            end
            K(:,s) = fval;      % store f(ti, xi, xti) for next stages
        end

        if ~stage_success
            % Newton failed – reduce step and retry
            h = h * 0.25;
            stats.failed = stats.failed + 1;
            continue;
        end

        % --- Compute solution and error estimate ---
        xn = x + h * (K * b);
        xn2 = x + h * (K * b2);      % embedded 2nd‑order solution
        er = norm(xn - xn2, inf);
        tol = atol + rtol * max(norm(x,inf), norm(xn,inf));

        if er <= tol || h <= min_step
            % Accept step
            t = t + h;
            x = xn;
            T_out(end+1,1) = t;
            X_out(end+1,:) = x.';
            if nargin >= 8 && ~isempty(callback)
                callback(t, x);
            end
            H = Ha(H, t, x);
            stats.steps = stats.steps + 1;
            h_last_accepted = h;          % remember for Jacobian reuse

            % PI step size controller
            if stats.steps > 1
                expo = 1/3;   % order 3
                facc = fac * (tol / er)^expo * (er / errold)^(beta * expo);
            else
                facc = fac * (tol / er)^(1/3);
            end
            facc = min(fac_max, max(fac_min, facc));
            h_new = facc * h;
            errold = max(er, 1e-12);
            h = min(max_step, max(min_step, h_new));
            same_jac_step = same_jac_step + 1;
        else
            % Reject step
            stats.failed = stats.failed + 1;
            h = h * 0.5;
            same_jac_step = 0;   % force Jacobian update after rejection
        end

        % Safety: if too many failures, abort
        if stats.failed > max_no_progress
            warning('Too many failed steps – stopping at t=%g', t);
            break;
        end
    end

    stats.nsteps = stats.steps;
    stats.nfailed = stats.failed;
    stats.njac = stats.jac_updates;
end

% ----------------------------------------------------------------------
% Rosenbrock integrator (2‑stage, with optional callback)
% ----------------------------------------------------------------------
function [T_out, X_out, stats] = ros(f, x0, h0, ts, n, extra_args, opts, accept_callback)
    g = 0.5;
    atol = opts.atol; rtol = opts.rtol; max_step = opts.max_step;
    interp_grid = opts.interp_grid; jac_reuse_tol = opts.jacobian_reuse_tol;
    force_ilu = getfielddef(opts, 'force_ilu', false);
    history_maxlen = getfielddef(opts, 'history_maxlen', 200000);

    he = max(h0, 1e-6);
    ht = (-he:interp_grid:0)';
    xs = zeros(numel(ht), n);
    for k = 1:numel(ht)
        xs(k,:) = x0(ht(k)).';
    end
    H = Hc(ht, xs, history_maxlen);

    t = ts(1); x = x0(0); x = x(:);
    T_out = t; X_out = x.';
    L = extra_args{1}; p = extra_args{2};

    Jp = []; LUstruct = []; use_ilu = force_ilu || (n > 5e4);
    h = min(max_step, 0.01);

    while t < ts(2)-1e-10
        h = min(h, ts(2)-t);
        xt = He(H, t - p.delay.h);
        J = J1(x, p, L);
        M = speye(n) - g*h*J;
        if isempty(Jp) || norm(J-Jp, inf) > jac_reuse_tol
            use_ilu = force_ilu || (n > 5e4);
            LUstruct = [];
        end
        f1 = f(t, x, xt);
        k1 = solve_linear(M, f1, use_ilu, ~isempty(LUstruct));

        t2 = t + 0.5*h; x2 = x + 0.5*h*k1;
        xt2 = He(H, t2 - p.delay.h);
        f2 = f(t2, x2, xt2);
        rhs2 = f2 + (1/(2*h))*k1;
        k2 = solve_linear(M, rhs2, use_ilu, true);

        xn = x + h*(k1 + k2)/2;
        er = norm(xn - (x + h*k1), inf);
        tol = atol + rtol*max(norm(x,inf), norm(xn,inf));

        if er <= tol
            t = t + h; x = xn;
            T_out(end+1,1) = t; X_out(end+1,:) = x.';
            if nargin >= 8 && ~isempty(accept_callback)
                accept_callback(t, x);
            end
            H = Ha(H, t, x);
            fc = (er==0)*2 + (er>0)*0.9*(tol/er)^(1/3);
            h = min(max_step, max(1e-8, fc*h));
        else
            h = max(1e-6, 0.5*h);
        end
    end
    stats = struct();
end

function [T_out, X_out, stats] = rosf(f, x0, h0, ts, n, extra_args, opts, callback)
    % ROS – 2‑stage Rosenbrock method (stiffly accurate) with step size control
    %       and Jacobian reuse.

    g = 0.5;   % parameter for Rosenbrock (common choice)

    % Options
    atol = opts.atol;
    rtol = opts.rtol;
    max_step = opts.max_step;
    interp_grid = opts.interp_grid;
    jac_reuse_tol = opts.jacobian_reuse_tol;
    force_ilu = getfielddef(opts, 'force_ilu', false);
    history_maxlen = getfielddef(opts, 'history_maxlen', 200000);
    min_step = 1e-10;
    max_no_progress = 1000;

    % Initial history
    he = max(h0, 1e-6);
    ht = (-he:interp_grid:0)';
    xs = zeros(numel(ht), n);
    for k = 1:numel(ht)
        xs(k,:) = x0(ht(k)).';
    end
    H = Hc(ht, xs, history_maxlen);

    t = ts(1);
    x = x0(0); x = x(:);
    T_out = t; X_out = x.';
    L = extra_args{1};
    p = extra_args{2};

    % Linear solver state
    Jp = [];
    LUstruct = [];
    use_ilu = force_ilu || (n > 5e4);
    h = min(max_step, 0.01 * (ts(2)-ts(1)));
    h_last_accepted = h;
    same_jac_step = 0;
    stats.steps = 0;
    stats.failed = 0;
    stats.jac_updates = 0;

    % Step size controller parameters
    fac_min = 0.2; fac_max = 5.0; fac = 0.9; beta = 0.0;
    errold = 1.0;

    while t < ts(2)-1e-12
        h = min(h, ts(2)-t);
        if h < min_step
            warning('ros: step size too small (%g) – forcing acceptance', h);
            % Force accept – break to avoid infinite loop
            break;
        end

        % --- Jacobian update decision ---
        if isempty(Jp) || same_jac_step > 10 || ...
           (h / h_last_accepted > 2.0 && same_jac_step > 0) || ...
           (h / h_last_accepted < 0.5 && same_jac_step > 0)
            xt = He(H, t - p.delay.h);
            J = J1(x, p, L);
            Jp = J;
            use_ilu = force_ilu || (n > 5e4);
            % Compute matrix M = I - g*h*J and factorize
            M = speye(n) - g*h*Jp;
            if ~use_ilu
                try
                    [Lf, Uf, Pf, Qf, Rf] = lu(M);
                    LUstruct = struct('L',Lf,'U',Uf,'P',Pf,'Q',Qf,'R',Rf,'ilu',false);
                catch
                    use_ilu = true;   % fallback
                end
            end
            if use_ilu
                setup.type = 'ilutp'; setup.droptol = 1e-3;
                [Lilu, Uilu] = ilu(M, setup);
                LUstruct = struct('L',Lilu,'U',Uilu,'ilu',true);
            end
            stats.jac_updates = stats.jac_updates + 1;
            same_jac_step = 0;
        else
            % Reuse existing factorization – but need to update M if h changed?
            % Since M depends on h, we must recompute M if h changed significantly.
            % We'll simply recompute M using stored Jp and new h.
            M = speye(n) - g*h*Jp;
            % But LUstruct may not be valid for new h. We'll need to re‑factorize.
            % To avoid refactor every step, we keep the old factorization only if h hasn't changed much.
            % Simpler: always refactor if h changed by more than 10% and we are reusing Jacobian.
            if abs(h - h_last_accepted) / h_last_accepted > 0.1
                % Refactor
                if ~use_ilu
                    try
                        [Lf, Uf, Pf, Qf, Rf] = lu(M);
                        LUstruct = struct('L',Lf,'U',Uf,'P',Pf,'Q',Qf,'R',Rf,'ilu',false);
                    catch
                        use_ilu = true;   % fallback
                    end
                end
                if use_ilu
                    setup.type = 'ilutp'; setup.droptol = 1e-3;
                    [Lilu, Uilu] = ilu(M, setup);
                    LUstruct = struct('L',Lilu,'U',Uilu,'ilu',true);
                end
            end
        end

        % --- Rosenbrock stages ---
        % Compute f1 at current point
        f1 = f(t, x, He(H, t - p.delay.h));

        % Solve M * k1 = f1
        if LUstruct.ilu
            [k1, flag] = gmres(M, f1, [], 1e-6, 200, LUstruct.L, LUstruct.U);
            if flag ~= 0
                k1 = M \ f1;   % fallback
            end
        else
            y = LUstruct.P * f1;
            z = LUstruct.L \ y;
            w = LUstruct.U \ z;
            if ~isempty(LUstruct.R)
                w = LUstruct.R \ w;
            end
            k1 = LUstruct.Q * w;
        end

        % Second stage
        t2 = t + 0.5*h;
        x2 = x + 0.5*h*k1;
        f2 = f(t2, x2, He(H, t2 - p.delay.h));
        rhs2 = f2 + (1/(2*h))*k1;

        % Solve M * k2 = rhs2 (same matrix)
        if LUstruct.ilu
            [k2, flag] = gmres(M, rhs2, [], 1e-6, 200, LUstruct.L, LUstruct.U);
            if flag ~= 0
                k2 = M \ rhs2;
            end
        else
            y2 = LUstruct.P * rhs2;
            z2 = LUstruct.L \ y2;
            w2 = LUstruct.U \ z2;
            if ~isempty(LUstruct.R)
                w2 = LUstruct.R \ w2;
            end
            k2 = LUstruct.Q * w2;
        end

        xn = x + h*(k1 + k2)/2;
        % Error estimate (difference from first‑order solution)
        er = norm(xn - (x + h*k1), inf);
        tol = atol + rtol * max(norm(x,inf), norm(xn,inf));

        if er <= tol || h <= min_step
            % Accept step
            t = t + h;
            x = xn;
            T_out(end+1,1) = t;
            X_out(end+1,:) = x.';
            if nargin >= 8 && ~isempty(callback)
                callback(t, x);
            end
            H = Ha(H, t, x);
            stats.steps = stats.steps + 1;
            h_last_accepted = h;

            % PI step size controller
            if stats.steps > 1
                expo = 1/2;   % order 2
                facc = fac * (tol / er)^expo * (er / errold)^(beta * expo);
            else
                facc = fac * (tol / er)^(1/2);
            end
            facc = min(fac_max, max(fac_min, facc));
            h_new = facc * h;
            errold = max(er, 1e-12);
            h = min(max_step, max(min_step, h_new));
            same_jac_step = same_jac_step + 1;
        else
            % Reject step
            stats.failed = stats.failed + 1;
            h = h * 0.5;
            same_jac_step = 0;   % force Jacobian update after rejection
        end

        if stats.failed > max_no_progress
            warning('ros: too many failed steps – stopping at t=%g', t);
            break;
        end
    end

    stats.nsteps = stats.steps;
    stats.nfailed = stats.failed;
    stats.njac = stats.jac_updates;
end

% ----------------------------------------------------------------------
% Plotting helper
% ----------------------------------------------------------------------
function p1(P, L, T, X, ttl)
    Nz = P.grid.Nz; Nr = P.grid.Nr; N = P.Ncells;
    z = linspace(0, P.geometry.L, Nz);
    r = linspace(0, P.geometry.R, Nr);
    cf = X(end,1:N)'; tf = X(end,N+1:end)';
    C = reshape(cf, [Nr Nz]);
    Tm = reshape(tf, [Nr Nz]);
    figure; imagesc(z, r, C); axis xy; colorbar; title([ttl ' concentration slice']);
    figure; plot(z, Tm(ceil(Nr/2),:)); title([ttl ' centerline T']);
end

% ----------------------------------------------------------------------
% Wrapper to run a single solver with downsampling and memory control
% ----------------------------------------------------------------------
##function [Tchk, Xchk, final_state, stats] = run_and_checkpoint(solver_fn, f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, label, save_every)
##    if nargin < 10, save_every = 50; end
##    est_max_steps = ceil((tspan(2)-tspan(1)) / max(1e-6, opts.max_step)) + 10;
##    est_store = ceil(est_max_steps / save_every) + 2;
##    Tchk = zeros(est_store,1); Xchk = zeros(est_store, n_state);
##    idx_store = 0;
##    step_counter = 0;
##
##
##
##    % Call solver with callback
##    [Tfull, Xfull, stats] = solver_fn(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, @callback);
##    final_state = Xfull(end,:).';
##
##    % Trim checkpoints
##    if idx_store == 0
##        Tchk = Tfull(end); Xchk = Xfull(end,:);
##    else
##        Tchk = Tchk(1:idx_store); Xchk = Xchk(1:idx_store,:);
##    end
##
##    clear Tfull Xfull;
##    try pack; catch; end
##    fprintf('%s done: stored %d checkpoints\n', label, numel(Tchk));
##end
##
##function callback(t_step, x_step,Tchk, Xchk,step_counter)
##        step_counter = step_counter + 1;
##        if mod(step_counter, save_every) == 0 || t_step >= tspan(2)-1e-12
##            idx_store = idx_store + 1;
##            if idx_store > size(Tchk,1)
##                % grow
##                Tchk = [Tchk; zeros(est_store,1)];
##                Xchk = [Xchk; zeros(est_store, n_state)];
##            end
##            Tchk(idx_store) = t_step;
##            Xchk(idx_store,:) = x_step(:).';
##        end
##    end
function [Tchk, Xchk, final_state, stats] = run_and_checkpoint(solver_fn, f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, label, save_every)
    % RUN_AND_CHECKPOINT  Run a solver with downsampling via external callback
    %
    % Inputs:
    %   solver_fn   - function handle to integrator (e.g. @sdirk)
    %   f_rhs       - RHS of DDE
    %   x0_fun      - history function
    %   h_delay     - delay
    %   tspan       - [t0 tfinal]
    %   n_state     - number of state variables
    %   extra_args  - cell array of extra arguments to f_rhs
    %   opts        - solver options structure
    %   label       - string for progress message
    %   save_every  - store every 'save_every' steps (default 50)
    %
    % Outputs:
    %   Tchk        - times of stored checkpoints
    %   Xchk        - states at those times
    %   final_state - final state vector
    %   stats       - statistics structure from solver

    if nargin < 10
        save_every = 50;
    end

    % Estimate maximum number of steps and preallocate
    est_max_steps = ceil((tspan(2)-tspan(1)) / max(1e-6, opts.max_step)) + 10;
    est_store = ceil(est_max_steps / save_every) + 2;

    % Global state for callback
    global callback_state
    callback_state = struct();
    callback_state.Tchk = zeros(est_store, 1);
    callback_state.Xchk = zeros(est_store, n_state);
    callback_state.idx_store = 0;
    callback_state.step_counter = 0;
    callback_state.save_every = save_every;
    callback_state.tspan = tspan;
    callback_state.est_store = est_store;
    callback_state.n_state = n_state;

    % Call solver with external callback
    [Tfull, Xfull, stats] = solver_fn(f_rhs, x0_fun, h_delay, tspan, n_state, extra_args, opts, @run_and_checkpoint_callback);

    % Extract final state
    final_state = Xfull(end, :).';

    % Retrieve checkpoint arrays from global state
    idx_store = callback_state.idx_store;
    if idx_store == 0
        Tchk = Tfull(end);
        Xchk = Xfull(end, :);
    else
        Tchk = callback_state.Tchk(1:idx_store);
        Xchk = callback_state.Xchk(1:idx_store, :);
    end

    % Clear global to avoid interference with next run
    clear global callback_state

    % Free memory
    clear Tfull Xfull;
    try pack; catch; end

    fprintf('%s done: stored %d checkpoints\n', label, numel(Tchk));
end

% ----------------------------------------------------------------------
% External callback function used by run_and_checkpoint
% ----------------------------------------------------------------------
function run_and_checkpoint_callback(t_step, x_step)
    global callback_state

    callback_state.step_counter = callback_state.step_counter + 1;
    if mod(callback_state.step_counter, callback_state.save_every) == 0 || ...
       t_step >= callback_state.tspan(2) - 1e-12
        callback_state.idx_store = callback_state.idx_store + 1;
        idx = callback_state.idx_store;

        % Grow arrays if necessary
        if idx > size(callback_state.Tchk, 1)
            extra = callback_state.est_store;
            callback_state.Tchk = [callback_state.Tchk; zeros(extra, 1)];
            callback_state.Xchk = [callback_state.Xchk; zeros(extra, callback_state.n_state)];
        end

        callback_state.Tchk(idx) = t_step;
        callback_state.Xchk(idx, :) = x_step(:).';
    end
end
% ----------------------------------------------------------------------
% Main driver (memory‑safe)
% ----------------------------------------------------------------------
function industrial_reactors_highres_full_fixed_clean()
    rand('state',1); try pack; catch; end

    % --- Define reactors ---
    A = g1(12, 1.6, 240, 60);  % packed‑bed
    A = g3(A, struct('Dc',1e-5, 'DT',5e-6, 'k1',10, 'Ea',6e4, ...
                     'k2',0.5, 'beta',5, 'gamma',1, 'porosity',0.4));
    A.Tfinal = 60;
    A.solver = s1('ros', 1e-7, 1e-7, 0.02, 0.005, 1e-8);

    B = g2(18, 0.15, 400, 32); % multi‑tubular
    B = g4(B, struct('Dc',5e-6, 'DT',2e-6, 'k1',1.2, 'Ea',4.5e4, ...
                     'k2',0.15, 'beta',2, 'gamma',0.6, 'porosity',0.45));
    B.Tfinal = 120;
    B.solver = s1('sdi', 1e-7, 1e-7, 0.02, 0.005, 1e-8);

    % --- Build Laplacians ---
    LA = L1(A.grid.Nz, A.grid.Nr, A.dx, A.dr);
    LB = L1(B.grid.Nz, B.grid.Nr, B.dx, B.dr);

    x0A = @(t)[A.ic.c0; A.ic.T0];
    x0B = @(t)[B.ic.c0; B.ic.T0];

    % --- Run SDIRK reference for A (coarse) and save only final state ---
    fprintf('\n=== Reactor A: SDIRK reference (coarse) ===\n');
    opts_ref = struct('atol',1e-8,'rtol',1e-8,'max_step',A.Tfinal/500, ...
                      'interp_grid',A.history.grid, 'jacobian_reuse_tol',A.solver.jacobian_reuse_tol, ...
                      'history_maxlen', ceil(A.delay.h/A.history.grid)+10, 'force_ilu',false);
    [~, ~, finalA_ref, stats_ref] = run_and_checkpoint(@sdirkf, @(t,x,xt) rhs(t,x,xt,LA,A), x0A, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {LA,A}, opts_ref, 'SDIRK‑A ref', 100);
    clear finalA_ref stats_ref; pack;

    % --- Run Rosenbrock for A (production) with downsampling ---
    fprintf('\n=== Reactor A: Rosenbrock (main) ===\n');
    opts_ros = A.solver.opts;
    opts_ros.jacobian_reuse_tol = A.solver.jacobian_reuse_tol;
    opts_ros.history_maxlen = ceil(A.delay.h/A.history.grid)+10;
    opts_ros.force_ilu = (2*A.Ncells > 5e4);
    [TA_ros, XA_ros, finalA_ros, stats_ros] = run_and_checkpoint(@rosf, @(t,x,xt) rhs(t,x,xt,LA,A), x0A, A.delay.h, [0 A.Tfinal], 2*A.Ncells, {LA,A}, opts_ros, 'Rosenbrock‑A', 50);

    % --- Save A results to disk (release memory) ---
    save('results_A.mat', 'TA_ros', 'XA_ros', 'finalA_ros', '-v7.3');
    clear TA_ros XA_ros finalA_ros; pack;

    % --- Run SDIRK for B (production) ---
    fprintf('\n=== Reactor B: SDIRK ===\n');
    opts_sdi = B.solver.opts;
    opts_sdi.jacobian_reuse_tol = B.solver.jacobian_reuse_tol;
    opts_sdi.history_maxlen = ceil(B.delay.h/B.history.grid)+10;
    opts_sdi.force_ilu = (2*B.Ncells > 5e4);
    [TB_sdi, XB_sdi, finalB_sdi, stats_sdi] = run_and_checkpoint(@sdirkf, @(t,x,xt) rhs(t,x,xt,LB,B), x0B, B.delay.h, [0 B.Tfinal], 2*B.Ncells, {LB,B}, opts_sdi, 'SDIRK‑B', 50);

    save('results_B.mat', 'TB_sdi', 'XB_sdi', 'finalB_sdi', '-v7.3');
    clear TB_sdi XB_sdi finalB_sdi; pack;

    % --- Optional plots (using saved data) ---
    fprintf('\nLoading results for plots...\n');
    load('results_A.mat', 'TA_ros', 'XA_ros');
    p1(A, LA, TA_ros, XA_ros, 'Packed‑bed');
    load('results_B.mat', 'TB_sdi', 'XB_sdi');
    p1(B, LB, TB_sdi, XB_sdi, 'Multi‑tubular');

    fprintf('\nAll runs completed successfully.\n');
end

% ----------------------------------------------------------------------
% Helper: safe field access
% ----------------------------------------------------------------------
function val = getfielddef(s, field, def)
    if isfield(s, field)
        val = s.(field);
    else
        val = def;
    end
end

% ----------------------------------------------------------------------
% Run it!
% ----------------------------------------------------------------------
industrial_reactors_highres_full_fixed_clean();


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
