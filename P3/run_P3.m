%% run_P3.m  —  P3 (a)(c)(e)(f) + 보너스 일괄 실행
% 시간: 3 s 실행 약 40회 + 20 s 실행 6회 (보너스). 보너스는 맨 아래 do_bonus로 끌 수 있음.

clear; clc; close all;
do_bonus = false;

%% 0. 모델 지정 -------------------------------------------------------------
mdl = '';   % 비워 두면 열려 있는 모델 자동 선택. 직접 쓰려면 예: 'P3_model'
if isempty(mdl)
    cand = find_system('SearchDepth',0,'Type','block_diagram','BlockDiagramType','model');
    if numel(cand) ~= 1
        error('열려 있는 모델이 %d개입니다. mdl = ''P3_model''; 처럼 직접 적어 주세요.', numel(cand));
    end
    mdl = cand{1};
end
load_system(mdl);
mws = get_param(mdl, 'ModelWorkspace');
fprintf('모델: %s\n', mdl);
diary('P3_results.txt');

lam = 40; Jm = 0.01; Ts = 1e-3; B = 0.10; Bm = 0.03;
base = struct('beta',0.37, 'phi',0.05, 'rho',0, 'adapt',0, 'dob',0, ...
              'dz',0, 'vd',0, 'sig',0);

%% 1. (c) 경계층 비교 --------------------------------------------------------
fprintf('실행 중...\n');
phis = [0 0.02 0.05 0.2];
C = cell(1,4);
for i = 1:4
    C{i} = runcase(mdl, mws, withp(base,'phi',phis(i)), 3);
end
fprintf('  (c) done\n');

%% 2. (e) ASMC ---------------------------------------------------------------
E1 = runcase(mdl, mws, withp(base,'beta',0.2), 3);                       % SMC 작은 beta
E2 = C{3};                                                               % SMC beta=0.37, phi=0.05
rhos = [0 50 500];
EA = cell(1,3);
for i = 1:3
    EA{i} = runcase(mdl, mws, withp(base,'beta',0.2,'adapt',1,'rho',rhos(i)), 3);
end
fprintf('  (e) done\n');

%% 3. (f) DOB beta sweep -----------------------------------------------------
betas = 0:0.05:0.7;
EF = zeros(2, numel(betas));
for d = 0:1
    for i = 1:numel(betas)
        o = runcase(mdl, mws, withp(base,'beta',betas(i),'dob',d), 3);
        EF(d+1,i) = emax(o, 1, 3);
    end
end
fprintf('  (f) sweep done\n');
bmin = nan(1,2);
for d = 1:2
    j = find(EF(d,:) < 0.5e-3, 1, 'first');
    if ~isempty(j), bmin(d) = betas(j); end
end
F_sat = cell(1,2);  F_sgn = cell(1,2);
for d = 1:2
    F_sat{d} = runcase(mdl, mws, withp(base,'beta',bmin(d),'dob',d-1), 3);
    F_sgn{d} = runcase(mdl, mws, withp(base,'beta',bmin(d),'dob',d-1,'phi',0), 3);
end
fprintf('  (f) done\n');

%% 4. 보너스: 잡음 + 위치 차분 속도, 데드존 --------------------------------
if do_bonus
    dzs = [0 0.1 0.2 0.3 0.5];
    pb  = withp(base,'beta',0.2,'adapt',1,'rho',0);
    X0  = runcase(mdl, mws, pb, 20);                                     % 잡음 없음 기준
    X   = cell(1, numel(dzs));
    for i = 1:numel(dzs)
        X{i} = runcase(mdl, mws, withp(pb,'sig',1e-4,'vd',1,'dz',dzs(i)), 20);
        fprintf('  bonus dz=%.1f done\n', dzs(i));
    end
end

runcase_reset(mws, base);

%% 5. 표 출력 ----------------------------------------------------------------
fprintf('\n===== (c) 경계층 (R2, beta = 0.37) =====\n');
fprintf('규칙 beta*Ts/Jm = %.3f rad/s\n', 0.37*Ts/Jm);
fprintf('%-8s | max|e| t>1 [mrad] | 채터링 지수 | 규칙 | 2phi/lam [mrad] | 2phi/(e*lam) [mrad] | max|u|\n','제어');
nm = {'sgn','phi=0.02','phi=0.05','phi=0.2'};
for i = 1:4
    ph = phis(i);
    if ph == 0
        rule = '-'; b1 = NaN; b2 = NaN;
    else
        if ph >= 0.37*Ts/Jm, rule = '만족'; else, rule = '미달'; end
        b1 = 2*ph/lam*1e3; b2 = 2*ph/(exp(1)*lam)*1e3;
    end
    fprintf('%-8s |      %7.3f      |   %7.4f   | %-4s |     %6.2f      |       %6.2f        | %5.3f\n', ...
        nm{i}, 1e3*emax(C{i},1,3), chat(C{i}), rule, b1, b2, max(abs(sq(C{i}.u.Data))));
end

fprintf('\n===== (e) ASMC (R2, beta = 0.2, phi = 0.05, k = 1e-3) =====\n');
fprintf('%-14s | max|e| [0,0.5) | [0.5,1.5) | [1.5,3) [mrad] | Bhat(0.5) | Bhat(1.49) | Bhat(3) | t90 [s] | 채터링\n','실행');
R = {E1, E2, EA{1}, EA{2}, EA{3}};
rn = {'SMC b=0.2','SMC b=0.37','ASMC rho=0','ASMC rho=50','ASMC rho=500'};
for i = 1:5
    o = R{i};
    fprintf('%-14s |    %7.3f     |  %7.3f  |    %7.3f     |   %5.3f   |   %5.3f    |  %5.3f  | %6.3f  | %.4f\n', rn{i}, ...
        1e3*emax(o,0,0.5), 1e3*emax(o,0.5,1.5), 1e3*emax(o,1.5,3), ...
        Bat(o,0.5), Bat(o,1.49), Bat(o,3), t90(o, Bm, B), chat(o));
end
fprintf('Bhat 최댓값 (rho=0): %.3f at t = %.2f s\n', max(sq(EA{1}.Bhat.Data)), ...
    EA{1}.Bhat.Time(find(sq(EA{1}.Bhat.Data) == max(sq(EA{1}.Bhat.Data)), 1)));

fprintf('\n===== (f) DOB 유무, 최소 beta (phi = 0.05, 기준 0.5 mrad) =====\n');
fprintf('beta   | DOB 없음 max|e| [mrad] | DOB 있음 max|e| [mrad]\n');
for i = 1:numel(betas)
    fprintf('%4.2f   |        %8.3f         |        %8.3f\n', betas(i), 1e3*EF(1,i), 1e3*EF(2,i));
end
fprintf('최소 beta: DOB 없음 %.2f, DOB 있음 %.2f\n', bmin(1), bmin(2));
fprintf('채터링 지수 (sat, phi=0.05): DOB 없음 %.4f, DOB 있음 %.4f\n', chat(F_sat{1}), chat(F_sat{2}));
fprintf('채터링 지수 (sgn, 같은 beta): DOB 없음 %.4f, DOB 있음 %.4f\n', chat(F_sgn{1}), chat(F_sgn{2}));

if do_bonus
    fprintf('\n===== 보너스: x 잡음 1e-4 rad RMS, 속도 = 위치 차분, rho = 0 =====\n');
    s_noise = std(sq(X{1}.s.Data(X{1}.s.Time > 2)));
    fprintf('s 잡음 크기 (std, dz=0, t>2): %.3f rad/s\n', s_noise);
    fprintf('%-10s | Bhat(2) | Bhat(10) | Bhat(20) | 참 e RMS t>2 [mrad]\n','데드존');
    fprintf('%-10s |  %5.3f  |  %5.3f   |  %5.3f   |  %6.3f\n', '잡음없음', Bat(X0,2), Bat(X0,10), Bat(X0,20), 1e3*etrue_rms(X0));
    for i = 1:numel(dzs)
        fprintf('dz = %-5.1f |  %5.3f  |  %5.3f   |  %5.3f   |  %6.3f\n', dzs(i), ...
            Bat(X{i},2), Bat(X{i},10), Bat(X{i},20), 1e3*etrue_rms(X{i}));
    end
end
diary off;

%% 6. 그림 -------------------------------------------------------------------
% (a) 면 위의 오차 (해석식)
fa = figure('Name','P3 (a) e(t) on surface','Position',[40 40 700 300]);
tt = 0:1e-4:0.2;
plot(tt*1e3, (1 - lam*tt).*exp(-lam*tt), 'LineWidth', 1.2); hold on; grid on
xline(1e3/lam, '--', 't = 1/\lambda');  xline(2e3/lam, ':', 't = 2/\lambda');
yline(0, 'k');
xlabel('t [ms]'); ylabel('e(t) / e_0'); title('s \equiv 0, \lambda = 40 rad/s')
savepng(fa, 'fig_P3a.png');

% (c) 경계층
cl = lines(4);
fc = figure('Name','P3 (c) boundary layer','Position',[60 60 900 850]);
subplot(3,1,1); hold on; grid on
for i = 1:4, plot(C{i}.e.Time, 1e3*sq(C{i}.e.Data), 'Color', cl(i,:)); end
ylabel('e [mrad]'); xlabel('t [s]'); legend(nm, 'Location','eastoutside')
subplot(3,1,2); hold on; grid on
for i = [1 2 3 4], plot(C{i}.u.Time, sq(C{i}.u.Data), 'Color', cl(i,:)); end
ylabel('u [N\cdotm]'); xlabel('t [s]'); legend(nm, 'Location','eastoutside')
subplot(3,1,3); hold on; grid on
for i = 1:4, stairs(C{i}.u.Time, sq(C{i}.u.Data), 'Color', cl(i,:)); end
xlim([1.00 1.03]); ylabel('u [N\cdotm]'); xlabel('t [s]'); title('u 확대 (30 ms)')
legend(nm, 'Location','eastoutside')
savepng(fc, 'fig_P3c.png');

% (e) ASMC
ce = [0.5 0.5 0.5; 0 0 0; 0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188];
fe = figure('Name','P3 (e) ASMC','Position',[80 80 900 850]);
subplot(3,1,1)
for i = 1:5
    semilogy(R{i}.e.Time, 1e3*abs(sq(R{i}.e.Data)), 'Color', ce(i,:)); hold on
end
grid on; ylim([1e-3 30]); ylabel('|e| [mrad]'); xlabel('t [s]'); legend(rn, 'Location','eastoutside')
subplot(3,1,2); hold on; grid on
for i = 1:5, plot(R{i}.u.Time, sq(R{i}.u.Data), 'Color', ce(i,:)); end
ylabel('u [N\cdotm]'); xlabel('t [s]'); legend(rn, 'Location','eastoutside')
subplot(3,1,3); hold on; grid on
for i = 3:5, plot(R{i}.Bhat.Time, sq(R{i}.Bhat.Data), 'Color', ce(i,:), 'LineWidth', 1.2); end
yline(B, 'k--', 'B'); yline(Bm, 'k:', 'B_m');
ylabel('B hat [N\cdotm\cdots]'); xlabel('t [s]'); legend(rn(3:5), 'Location','eastoutside')
savepng(fe, 'fig_P3e.png');

% (f) beta sweep
ff = figure('Name','P3 (f) DOB beta sweep','Position',[100 100 700 350]);
semilogy(betas, 1e3*EF(1,:), 'o-', betas, 1e3*EF(2,:), 's-'); hold on; grid on
yline(0.5, 'r--', '0.5 mrad');
xlabel('\beta [N\cdotm]'); ylabel('max_{t>1}|e| [mrad]'); legend('DOB 없음','DOB 있음 (\omega_c = 100)')
savepng(ff, 'fig_P3f.png');

% 보너스
if do_bonus
    fb = figure('Name','P3 bonus drift','Position',[120 120 800 350]);
    plot(X0.Bhat.Time, sq(X0.Bhat.Data), 'k', 'LineWidth', 1.2); hold on; grid on
    for i = 1:numel(dzs), plot(X{i}.Bhat.Time, sq(X{i}.Bhat.Data)); end
    yline(B, 'k--', 'B');
    xlabel('t [s]'); ylabel('B hat [N\cdotm\cdots]')
    legend(['잡음 없음', arrayfun(@(z) sprintf('\\phi_a = %.1f', z), dzs, 'UniformOutput', false)], ...
        'Location','eastoutside')
    savepng(fb, 'fig_P3bonus.png');
end

if do_bonus
    save('P3_results.mat','C','E1','E2','EA','EF','betas','bmin','F_sat','F_sgn','X0','X','dzs','-v7.3');
else
    save('P3_results.mat','C','E1','E2','EA','EF','betas','bmin','F_sat','F_sgn','-v7.3');
end
fprintf('\n완료: P3_results.txt, P3_results.mat, fig_P3a/c/e/f(/bonus).png 저장\n');


%% ======================= 함수 =======================
function o = runcase(mdl, mws, p, stopT)
    P = struct('Jm',0.01,'Bm',0.03,'lam',40,'Ts',1e-3,'umax',2, ...
               'beta',p.beta,'phi',p.phi,'k',1e-3,'rho',p.rho,'wc',100, ...
               'adapt_on',p.adapt,'dob_on',p.dob,'dz',p.dz,'vel_diff',p.vd);
    assignin(mws, 'P', P);
    assignin(mws, 'sigma_x', p.sig);
    o = sim(mdl, 'StopTime', num2str(stopT));
end
function runcase_reset(mws, p)
    P = struct('Jm',0.01,'Bm',0.03,'lam',40,'Ts',1e-3,'umax',2, ...
               'beta',p.beta,'phi',p.phi,'k',1e-3,'rho',p.rho,'wc',100, ...
               'adapt_on',p.adapt,'dob_on',p.dob,'dz',p.dz,'vel_diff',p.vd);
    assignin(mws, 'P', P);  assignin(mws, 'sigma_x', 0);
end
function p = withp(p, varargin)
    for i = 1:2:numel(varargin), p.(varargin{i}) = varargin{i+1}; end
end
function v = sq(x), v = squeeze(x); end
function m = emax(o, t1, t2)
    t = o.e.Time;  v = sq(o.e.Data);
    m = max(abs(v(t >= t1 & t < t2)));
end
function c = chat(o)
    c = mean(abs(diff(sq(o.u.Data))));
end
function b = Bat(o, tq)
    b = interp1(o.Bhat.Time, sq(o.Bhat.Data), tq, 'previous');
end
function t = t90(o, Bm, B)
    v = sq(o.Bhat.Data);  j = find(v >= Bm + 0.9*(B - Bm), 1, 'first');
    if isempty(j), t = NaN; else, t = o.Bhat.Time(j); end
end
function r = etrue_rms(o)
    t = o.x.Time;  e = 0.5*(1 - cos(2*pi*t)) - sq(o.x.Data);
    r = sqrt(mean(e(t > 2).^2));
end
function savepng(f, name)
    try, exportgraphics(f, name, 'Resolution', 200); catch, saveas(f, name); end
end