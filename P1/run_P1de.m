%% run_P1de.m  —  P1 (d)(e) 일괄 실행 스크립트

clear; clc; close all;

%% 0. 모델 지정 -------------------------------------------------------------
mdl = '';   % 비워 두면 열려 있는 모델을 자동으로 찾음. 직접 쓰려면 예: 'P1_model'

if isempty(mdl)
    cand = find_system('SearchDepth', 0, 'Type', 'block_diagram', ...
                       'BlockDiagramType', 'model');
    if numel(cand) ~= 1
        error('열려 있는 모델이 %d개입니다. 위의 mdl = ''모델이름''; 에 직접 적어 주세요.', numel(cand));
    end
    mdl = cand{1};
end
load_system(mdl);
mws = get_param(mdl, 'ModelWorkspace');
fprintf('모델: %s\n', mdl);

diary('P1_results.txt');   % 명령 창 출력을 파일로도 저장

%% 1. 케이스 정의 ------------------------------------------------------------
base = struct('mode',1, 'zeta',0.7, 'wa',200, 'sigma',0, ...
              'Bn',0.06, 'umax',6, 'd_on',1);
dg   = withp(base, 'Bn',0.10, 'umax',1e9);      % 진단 실행: 감쇠 오차·포화 없음

%% 2. 시뮬레이션 ------------------------------------------------------------
fprintf('실행 중...\n');
% 실제 조건 (d)
R1  = runcase(mdl, mws, base);                                  fprintf('  R1 done\n');
R2  = runcase(mdl, mws, withp(base,'zeta',0.3));                fprintf('  R2 done\n');
R3  = runcase(mdl, mws, withp(base,'mode',2,'wa',200));         fprintf('  R3 done\n');
% 진단 실행 (정착시간용: 외란 있음 / 없음 쌍)
D1  = runcase(mdl, mws, dg);
D1b = runcase(mdl, mws, withp(dg,'d_on',0));                    fprintf('  D1 done\n');
D2  = runcase(mdl, mws, withp(dg,'zeta',0.3));
D2b = runcase(mdl, mws, withp(dg,'zeta',0.3,'d_on',0));         fprintf('  D2 done\n');
D3  = runcase(mdl, mws, withp(dg,'mode',2,'wa',200));
D3b = runcase(mdl, mws, withp(dg,'mode',2,'wa',200,'d_on',0));  fprintf('  D3 done\n');
% 노이즈 (e)
R4  = runcase(mdl, mws, withp(base,'sigma',0.05));                       fprintf('  R4 done\n');
R5  = runcase(mdl, mws, withp(base,'mode',2,'wa',200,'sigma',0.05));     fprintf('  R5 done\n');
R6  = runcase(mdl, mws, withp(base,'mode',2,'wa',1000,'sigma',0.05));    fprintf('  R6 done\n');

runset(mws, base);   % 기본값으로 복구

%% 3. (d) 표 ----------------------------------------------------------------
names = {'conv zeta=0.7', 'conv zeta=0.3', 'TI wa=200'};
R  = {R1, R2, R3};
D  = {D1, D2, D3};
Db = {D1b, D2b, D3b};

fprintf('\n===== (d) max|s| [rad/s] 와 d_hat 정착시간 =====\n');
fprintf('%-14s | 기저[0.5,1) | t=1후 0.5s | 기저[1.5,2) | t=2후 0.5s | 정착[s] | 진단 max|s| t=1후\n', '케이스');
for i = 1:3
    fprintf('%-14s |   %7.3f   |  %7.3f   |   %7.3f   |  %7.3f   | %7.4f |  %7.4f\n', names{i}, ...
        maxs(R{i},0.5,1), maxs(R{i},1,1.5), maxs(R{i},1.5,2), maxs(R{i},2,2.5), ...
        settle(D{i}, Db{i}), maxs_step(D{i}, Db{i}));
end

fprintf('\n===== telescoped PI 형태 검증 (진단 실행 D3, 포화 없음) =====\n');
fprintf('max|d_hat - d_hat_tele| = %.3e\n', ...
    max(abs(sq(D3.d_hat.Data) - sq(D3.d_hat_tele.Data))));

%% 4. (e) 표 ----------------------------------------------------------------
E  = {R4, R5, R6};
en = {'conv zeta=0.7', 'TI wa=200', 'TI wa=1000'};
Kp = [0.02*30, 0.02*(30+200), 0.02*(30+1000)];

fprintf('\n===== (e) 노이즈 sigma = 0.05 rad/s =====\n');
fprintf('%-14s |  Kp   | sqrt2*Kp*sig | RMS(s)측정 | RMS(s)참 | RMS(du) | u간격[ms]\n', '케이스');
for i = 1:3
    o = E{i};
    s_meas = sq(o.s.Data);
    s_true = sq(o.x.Data) - 40*sin(2*pi*o.x.Time);
    du     = diff(sq(o.u.Data));
    fprintf('%-14s | %5.2f |    %6.3f    |  %7.4f   | %7.4f  | %7.4f |  %5.2f\n', en{i}, ...
        Kp(i), sqrt(2)*Kp(i)*0.05, rmsv(s_meas), rmsv(s_true), rmsv(du), ...
        1e3*median(diff(o.u.Time)));
end
fprintf('\n※ u간격이 1.00이 아니면 out.u가 Ts로 기록되지 않은 것 → RMS(du) 무효\n');

diary off;

%% 5. 그림 -------------------------------------------------------------------
c7 = [0 0.447 0.741];  c3 = [0.850 0.325 0.098];  cT = [0.466 0.674 0.188];

% (d)(iv) 실제 조건 비교: (c)와 같은 형식
f1 = figure('Name','(d) conventional vs transient-improved', 'Position',[100 100 900 800]);
subplot(3,1,1); hold on; grid on
plot(R1.s.Time, sq(R1.s.Data), 'Color', c7)
plot(R2.s.Time, sq(R2.s.Data), 'Color', c3)
plot(R3.s.Time, sq(R3.s.Data), 'Color', cT, 'LineWidth', 1.2)
ylabel('s [rad/s]'); legend('conv \zeta=0.7','conv \zeta=0.3','TI \omega_a=200','Location','best')

subplot(3,1,2); hold on; grid on
plot(R1.delta.Time, sq(R1.delta.Data), 'k')
plot(R1.d_hat.Time, sq(R1.d_hat.Data), 'Color', c7)
plot(R2.d_hat.Time, sq(R2.d_hat.Data), 'Color', c3)
plot(R3.d_hat.Time, sq(R3.d_hat.Data), 'Color', cT, 'LineWidth', 1.2)
ylabel('[N\cdotm]'); legend('\delta','conv \zeta=0.7','conv \zeta=0.3','TI \omega_a=200','Location','best')

subplot(3,1,3); hold on; grid on
plot(R1.u.Time, sq(R1.u.Data), 'Color', c7)
plot(R2.u.Time, sq(R2.u.Data), 'Color', c3)
plot(R3.u.Time, sq(R3.u.Data), 'Color', cT, 'LineWidth', 1.2)
ylabel('u [N\cdotm]'); xlabel('t [s]')
legend('conv \zeta=0.7','conv \zeta=0.3','TI \omega_a=200','Location','best')
savefig_png(f1, 'fig_d_compare.png');

% (d) 진단: t=1 스텝 성분만 (외란 있음 - 없음)
f2 = figure('Name','(d) 진단 스텝 응답', 'Position',[150 150 900 600]);
subplot(2,1,1); hold on; grid on
t = D1.d_hat.Time;
plot(t, 0.5*(t>=1), 'k')
plot(t, sq(D1.d_hat.Data)-sq(D1b.d_hat.Data), 'Color', c7)
plot(t, sq(D2.d_hat.Data)-sq(D2b.d_hat.Data), 'Color', c3)
plot(t, sq(D3.d_hat.Data)-sq(D3b.d_hat.Data), 'Color', cT, 'LineWidth', 1.2)
xlim([0.95 1.5]); ylabel('d\_hat 스텝 성분 [N\cdotm]')
legend('\delta (스텝)','conv \zeta=0.7','conv \zeta=0.3','TI \omega_a=200','Location','best')

subplot(2,1,2); hold on; grid on
plot(D1.s.Time, sq(D1.s.Data)-sq(D1b.s.Data), 'Color', c7)
plot(D2.s.Time, sq(D2.s.Data)-sq(D2b.s.Data), 'Color', c3)
plot(D3.s.Time, sq(D3.s.Data)-sq(D3b.s.Data), 'Color', cT, 'LineWidth', 1.2)
xlim([0.95 1.5]); ylabel('s 스텝 성분 [rad/s]'); xlabel('t [s]')
legend('conv \zeta=0.7','conv \zeta=0.3','TI \omega_a=200','Location','best')
savefig_png(f2, 'fig_d_diag.png');

save('P1_results.mat', 'R1','R2','R3','R4','R5','R6','D1','D1b','D2','D2b','D3','D3b');
fprintf('\n완료: P1_results.txt, P1_results.mat, fig_d_compare.png, fig_d_diag.png 저장\n');


%% ======================= 함수 =======================
function o = runcase(mdl, mws, p)
    runset(mws, p);
    o = sim(mdl);
end

function runset(mws, p)
    J = 0.02; lam = 30; Ts = 1e-3;
    assignin(mws, 'mode',    p.mode);
    assignin(mws, 'zeta',    p.zeta);
    assignin(mws, 'k_conv',  J*lam^2/(4*p.zeta^2));
    assignin(mws, 'wa',      p.wa);
    assignin(mws, 'ka',      p.wa*Ts);
    assignin(mws, 'sigma_n', p.sigma);
    assignin(mws, 'Bn',      p.Bn);
    assignin(mws, 'umax',    p.umax);
    assignin(mws, 'd_on',    p.d_on);
end

function p = withp(p, varargin)
    for i = 1:2:numel(varargin)
        p.(varargin{i}) = varargin{i+1};
    end
end

function v = sq(x)
    v = squeeze(x);
end

function r = rmsv(x)
    r = sqrt(mean(x.^2));
end

function m = maxs(o, t1, t2)
    t = o.s.Time;  v = sq(o.s.Data);
    m = max(abs(v(t >= t1 & t < t2)));
end

function m = maxs_step(o, ob)
    % 진단 실행에서 t=1 스텝 성분의 max|s| (기준선 제거)
    t = o.s.Time;  v = sq(o.s.Data) - sq(ob.s.Data);
    m = max(abs(v(t >= 1 & t < 1.5)));
end

function ts = settle(o, ob)
    % 외란 있음 - 없음 = 스텝 성분만. t=2 사인 시작 전 [1,2)에서 2% (0.01 N·m) 판정
    t = o.d_hat.Time;
    e = sq(o.d_hat.Data) - sq(ob.d_hat.Data) - 0.5;
    idx  = find(t >= 1 & t < 2);
    last = find(abs(e(idx)) > 0.01, 1, 'last');
    if isempty(last)
        ts = 0;
    else
        ts = t(idx(min(last+1, numel(idx)))) - 1;
    end
end

function savefig_png(f, name)
    try
        exportgraphics(f, name, 'Resolution', 200);
    catch
        saveas(f, name);
    end
end