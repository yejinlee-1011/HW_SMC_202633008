function [u, e, s, Bhat, dhat, u_cmd] = ctrl(xm, vs, xd, xd_dot, xd_ddot, P)
% ---- 출력 기본값 (모든 경로에서 정의되도록) ----
u = 0; e = 0; s = 0; Bhat = 0; dhat = 0; u_cmd = 0;

% ---- 상태 변수 ----
persistent ie Bh z xprev
if isempty(ie)
    ie = 0; Bh = P.Bm; z = 0; xprev = xm;
end

% ---- 속도: 측정값, 또는 위치 차분(보너스) ----
if P.vel_diff > 0.5
    vm = (xm - xprev)/P.Ts;
else
    vm = vs;
end
xprev = xm;

% ---- 슬라이딩 면 (e = xd - x) ----
e    = xd - xm;
de   = xd_dot - vm;
ie_t = ie + P.Ts*e;
s    = de + 2*P.lam*e + P.lam^2*ie_t;

% ---- 스위칭: phi = 0 이면 sgn, 아니면 sat ----
if P.phi == 0
    sw = sign(s);
else
    sw = min(max(s/P.phi, -1), 1);
end

% ---- 등가 제어 ----
if P.adapt_on > 0.5
    Bu = Bh;
else
    Bu = P.Bm;
end
ueq = P.Jm*xd_ddot + Bu*vm + P.Jm*(2*P.lam*de + P.lam^2*e);

% ---- DOB ----
dhat = P.dob_on*(z - P.Jm*P.wc*vm);

% ---- 제어 입력과 포화 ----
u_cmd = ueq + P.beta*sw + dhat;
u     = min(max(u_cmd, -P.umax), P.umax);

% ---- anti-windup: 포화 아닐 때만 적분 확정 ----
if abs(u_cmd) <= P.umax
    ie = ie_t;
end

% ---- 적응법칙 (누설 rho, 데드존 dz, clip [0, 0.5]) ----
if P.adapt_on > 0.5 && abs(s) >= P.dz
        Bh = Bh + P.Ts*P.k*(vm*sw/P.Jm - P.rho*(Bh - P.Bm));
    Bh = min(max(Bh, 0), 0.5);
end
Bhat = Bh;

% ---- DOB 상태 갱신 ----
z = z + P.Ts*P.wc*((u - P.Bm*vm + P.Jm*P.wc*vm) - z);