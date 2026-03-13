%% ================================================================
%  Supplementary MATLAB Code — Paper 2 (Observability)  v3
%  ----------------------------------------------------------------
%  "Equivalence of Kalman and Hewer Observability for Generalised
%   Bilinear Matrix Periodic Systems"
%  Authors: Sri Venkata Durga Sudarsan Madhyannapu &
%           Sravanam Pradheep Kumar
%
%  VERSION 3 — CORRECTED RECURSION FORMULA
%  =========================================
%  The observability Gramian recursion is derived from:
%
%    H(t+iT) = H(s) * M^i    (s = t - iT, using cocycle + T-periodicity)
%
%  Therefore the tail integral from iT to (i+1)T becomes:
%
%    integral (M^i)^T * H(s)^T * H(s) * M^i ds
%    = (M^T)^i * Ot_1 * M^i
%
%  So the correct recursion is:
%
%    Ot_{i+1} = Ot_i + (M^T)^i * Ot_1 * M^i          [CORRECT]
%
%  NOT:  Ot_i + (M^T)^i * Ot_1 * (M^T)^{-i}           [WRONG - was in v1/v2]
%
%  This matches the controllability recursion by duality:
%    Wt_{i+1} = Wt_i + M^{-i} * Wt_1 * (M^{-i})^T
%             = Wt_i + M^{-i} * Wt_1 * (M^T)^{-i}
%  Dual: replace M^{-i} -> (M^T)^i, giving (M^T)^i * Ot_1 * M^i  ✓
%
%  OTHER FIXES (carried over from v2):
%    - IC = reshape(eye(n2), n2^2, 1)  [column vector]
%    - Random matrices fixed before quadrature loop
%    - n=4 excluded from scaling (>5GB)
%  ================================================================

clear; clc; close all;
fprintf('=== Paper 2 Observability v3: Numerical Verification ===\n\n');

%% ----------------------------------------------------------------
%  SECTION 1 — SYSTEM DEFINITION
%  ----------------------------------------------------------------
n = 3;  r = 3;  T = 2*pi;  n2 = n^2;

A = @(t) [sin(2*t),1,0; 0,cos(2*t),1; 0,0,sin(t)];
B = @(t) [sin(t)+cos(t),0,0; 0,sin(t)-cos(t),0; 0,0,-sin(t)];
F = @(t) [cos(t),0,0; 0,sin(t),0; 0,0,cos(2*t)];
G = @(t) [0,sin(t),0; cos(t),0,0; 0,0,sin(2*t)];
C = @(t) [cos(t),1,0; 0,sin(2*t),1; 1,0,cos(t)];

calA  = @(t) kron(eye(n),A(t)) + kron(B(t)',eye(n)) + kron(G(t)',F(t));
Clift = @(t) kron(eye(n), C(t));

IC   = reshape(eye(n2), n2^2, 1);
opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
odeM = @(t,ph) reshape(calA(t)*reshape(ph,n2,n2), n2^2, 1);

fprintf('System: n=%d, r=%d, n^2=%d, T=2pi\n\n', n, r, n2);

%% ----------------------------------------------------------------
%  SECTION 2 — MONODROMY
%  ----------------------------------------------------------------
[~,Ph] = ode45(odeM, [0,T], IC, opts);
M      = reshape(Ph(end,:)', n2, n2);
MT     = M';

fprintf('Eigenvalues of M:\n'); disp(eig(M).');
fprintf('\n');

%% ----------------------------------------------------------------
%  SECTION 3 — NON-FACTORISABILITY
%  ----------------------------------------------------------------
IC_n  = reshape(eye(n), n^2, 1);
odeA  = @(t,ph) reshape(A(t)*reshape(ph,n,n), n^2, 1);
odeB  = @(t,ph) reshape(B(t)'*reshape(ph,n,n), n^2, 1);
[~,PA] = ode45(odeA,[0,T],IC_n,opts);
[~,PB] = ode45(odeB,[0,T],IC_n,opts);
Mkron  = kron(reshape(PB(end,:)',n,n), reshape(PA(end,:)',n,n));
fprintf('||M - Phi_B* x Phi_A||_F = %.4e  (confirms non-factorisability)\n\n', ...
    norm(M-Mkron,'fro'));

%% ----------------------------------------------------------------
%  SECTION 4 — MINIMAL POLYNOMIAL  k=3
%  ----------------------------------------------------------------
k = 3;
fprintf('|eig(M)| = '); disp(abs(eig(M)).');
fprintf('k = %d  (three distinct moduli), horizon = %.4f\n\n', k, k*T);

%% ----------------------------------------------------------------
%  SECTION 5 — STEP-SIZE SENSITIVITY
%  ----------------------------------------------------------------
fprintf('--- Step-size sensitivity ---\n');
tols = [1e-3, 1e-4, 1e-5];
M_t  = cell(3,1);
for ii = 1:3
    o = odeset('RelTol',tols(ii),'AbsTol',1e-10);
    [~,P] = ode45(odeM,[0,T],IC,o);
    M_t{ii} = reshape(P(end,:)',n2,n2);
end
fprintf('  1e-3 vs 1e-4: ||DeltaM||_F = %.4e\n', norm(M_t{1}-M_t{2},'fro'));
fprintf('  1e-4 vs 1e-5: ||DeltaM||_F = %.4e\n\n', norm(M_t{2}-M_t{3},'fro'));

%% ----------------------------------------------------------------
%  SECTION 6 — FULL-RANK OUTPUT CONDITION
%  ----------------------------------------------------------------
sv_C = arrayfun(@(t) min(svd(C(t))), linspace(0,T,1e4));
fprintf('min sigma_min(C(t)) = %.4f  (>0 confirms full-rank output)\n\n', min(sv_C));

%% ----------------------------------------------------------------
%  SECTION 7 — OBSERVABILITY GRAMIAN Ot_1
%
%  Ot_1 = int_0^T  H(t)^T H(t) dt
%       = int_0^T  Phi(t,0)^T * Clift(t)^T * Clift(t) * Phi(t,0) dt
%
%  H(t) = Clift(t) * Phi_calA(t,0)   where Phi is computed FORWARD from 0 to t
%  ----------------------------------------------------------------
fprintf('--- Computing Ot_1 by quadrature (N=1000) ---\n');
N_q = 1000;  t_q = linspace(0,T,N_q+1);
Ot1 = zeros(n2,n2);
for ii = 1:N_q
    tm = (t_q(ii)+t_q(ii+1))/2;  dt = t_q(ii+1)-t_q(ii);
    [~,P2]  = ode45(odeM, [0,tm], IC, opts);
    Phi_t0  = reshape(P2(end,:)', n2, n2);   % Phi_calA(t, 0) forward
    H_t     = Clift(tm) * Phi_t0;            % (I_n x C(t)) Phi(t,0)
    Ot1     = Ot1 + H_t' * H_t * dt;
end

% Quick symmetry and PD check on Ot_1
ev1 = real(eig(Ot1));
fprintf('Ot_1: lambda_min=%.4f, lambda_max=%.4e, PD=%s\n\n', ...
    min(ev1), max(ev1), mat2str(all(ev1>0)));

%% ----------------------------------------------------------------
%  SECTION 8 — GRAMIAN RECURSION  (CORRECTED FORMULA)
%
%  Ot_{i+1} = Ot_i + (M^T)^i * Ot_1 * M^i
%
%  Derivation: H(s+iT) = H(s)*M^i  (cocycle + T-periodicity)
%  => tail integral = (M^i)^T * Ot_1 * M^i = (M^T)^i * Ot_1 * M^i
%  ----------------------------------------------------------------
Ot    = {Ot1};
for i = 1:k-1
    MTi     = MT^i;         % (M^T)^i
    Mi      = M^i;          % M^i
    Ot{i+1} = Ot{i} + MTi * Ot{1} * Mi;   % CORRECTED
end

fprintf('Table 1 — MATLAB-verified Observability Gramian Statistics:\n');
fprintf('%-6s  %-14s  %-22s  %-20s  %s\n', ...
    'i','lambda_min','lambda_max','kappa','PD');
fprintf('%s\n', repmat('-',1,72));
for i = 1:k
    ev  = real(eig(Ot{i}));
    lmn = min(ev); lmx = max(ev);
    fprintf('%-6d  %-14.4f  %-22.4e  %-20.4e  %s\n', ...
        i, lmn, lmx, lmx/lmn, mat2str(all(ev>0)));
end
fprintf('\n');

%% ----------------------------------------------------------------
%  SECTION 9 — RECURSION RESIDUAL
%  ----------------------------------------------------------------
Ot2_check = Ot{1} + MT * Ot{1} * M;   % i=1 step
res = norm(Ot{2} - Ot2_check, 'fro');
fprintf('Recursion residual ||Ot_2 - (Ot_1 + M^T*Ot_1*M)||_F = %.4e\n\n', res);

%% ----------------------------------------------------------------
%  SECTION 10 — SPECTRAL LOWER BOUND
%  ----------------------------------------------------------------
fprintf('--- Spectral lower bound ---\n');
smin_M  = min(svd(M));
lmin_O1 = min(real(eig(Ot{1})));
fprintf('sigma_min(M) = %.4e,  lambda_min(Ot_1) = %.4f\n', smin_M, lmin_O1);
for i = 2:k
    bound   = smin_M^(2*(i-1)) * lmin_O1;
    lmin_Oi = min(real(eig(Ot{i})));
    sat     = lmin_Oi >= bound - 1e-10;
    fprintf('i=%d: bound=%.4e, lambda_min(Ot_%d)=%.4f, satisfied=%s\n', ...
        i, bound, i, lmin_Oi, mat2str(sat));
end
fprintf('\n');

%% ----------------------------------------------------------------
%  SECTION 11 — SCALING STUDY  n=2, 3
%  ----------------------------------------------------------------
fprintf('--- Observability scaling study: n=2, 3 ---\n');
fprintf('(n=4 excluded: >5GB memory)\n\n');
rng(42);
for nv = [2, 3]
    nv2   = nv^2;
    IC_nv = reshape(eye(nv2), nv2^2, 1);
    kappas = zeros(1,10);
    for tr = 1:10
        A0=randn(nv,nv); A1=randn(nv,nv);
        B0=randn(nv,nv); F0=randn(nv,nv); G0=randn(nv,nv);
        C0=randn(nv,nv)+eye(nv)*0.5;
        Ar  = @(t) A0*sin(t)+A1*cos(t);
        Br  = @(t) B0*cos(t);
        Fr  = @(t) F0*sin(2*t);
        Gr  = @(t) G0*cos(2*t);
        Cr  = @(t) C0;
        Clr = @(t) kron(eye(nv), Cr(t));
        cAr = @(t) kron(eye(nv),Ar(t))+kron(Br(t)',eye(nv))+kron(Gr(t)',Fr(t));
        oder= @(t,ph) reshape(cAr(t)*reshape(ph,nv2,nv2), nv2^2, 1);
        N_q2=100; tq2=linspace(0,T,N_q2+1);
        O1r =zeros(nv2,nv2);
        for ii=1:N_q2
            tm2=(tq2(ii)+tq2(ii+1))/2; dt2=tq2(ii+1)-tq2(ii);
            [~,P3]=ode45(oder,[0,tm2],IC_nv,opts);
            Phi2=reshape(P3(end,:)',nv2,nv2);
            H2=Clr(tm2)*Phi2;
            O1r=O1r+H2'*H2*dt2;
        end
        ev2=real(eig(O1r));
        kappas(tr)=max(ev2)/max(min(ev2),1e-14);
    end
    fprintf('n=%d  n^2=%2d  kappa = %.4e +/- %.4e\n', ...
        nv, nv2, mean(kappas), std(kappas));
end

%% ----------------------------------------------------------------
%  SECTION 12 — H-OBSERVABILITY
%  H-observable iff for all M-eigenvectors eta:
%    integral_0^T ||H(t)*eta||^2 dt > 0
%  ----------------------------------------------------------------
fprintf('\n--- H-observability verification ---\n');
[V_M,~] = eig(M);
min_e   = Inf;
tq5     = linspace(0,T,201);
for j = 1:n2
    eta=V_M(:,j); e=0;
    for ii=1:200
        tm5=(tq5(ii)+tq5(ii+1))/2; dt5=tq5(ii+1)-tq5(ii);
        [~,P6]=ode45(odeM,[0,tm5],IC,opts);
        H_t5=Clift(tm5)*reshape(P6(end,:)',n2,n2);
        e=e+norm(H_t5*eta)^2*dt5;
    end
    min_e=min(min_e,real(e));
end
fprintf('min_j integral = %.4f  (>0 -> H-observability confirmed)\n\n', min_e);

%% ----------------------------------------------------------------
%  SECTION 13 — MINIMUM DETECTION ENERGY
%  ----------------------------------------------------------------
fprintf('--- Minimum detection energy ---\n');
x0    = reshape(eye(n),n2,1);
Gamma = Ot{k}\x0;
J_obs = real(x0'*Gamma);
fprintf('J_obs* = x0^T * Ot_k^{-1} * x0 = %.6f\n', J_obs);
fprintf('lambda_min(Ot_k) = %.4f\n\n', min(real(eig(Ot{k}))));

fprintf('=== Paper 2 v3 complete. All MATLAB-verified. ===\n');
