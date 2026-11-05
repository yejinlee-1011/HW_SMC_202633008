%% run_P2.m  —  P2 (b)(c)(d)(e) + 보너스 일괄 실행

clear; clc; close all;

%% 0. 모델 지정 -------------------------------------------------------------
mdl = '';   % 비워 두면 열려 있는 모델 자동 선택. 직접 쓰려면 예: 'P2_model'
if isempty(mdl)
    cand = find_system('SearchDepth',0,'Type','block_diagram','BlockDiagramType','model');
    if numel(cand) ~= 1
        error('열려 있는 모델이 %d개입니다. mdl = ''P2_model''; 처럼 직접 적어 주세요.', numel(cand));
    end
    mdl = cand{1};
end
load_system(mdl);
mws = get_param(mdl, 'ModelWorkspace');
fprintf('모델: %s\n', mdl);
diary('P2_results.txt');

%% 1. 케이스 -----------------------------------------------------------------
base = struct('mode',1, 'zeta',0.7, 'wa',100, 'Jw',1.2, 'quant_on',0, 'act_on',0);

fprintf('실행 중...\n');
R1 = runcase(mdl, mws, base);                                   fprintf('  R1 (b) conv\n');
R2 = runcase(mdl, mws, withp(base,'mode',2));                   fprintf('  R2 (c) TI wa=100\n');
R3 = runcase(mdl, mws, withp(base,'mode',2,'Jw',1.44));         fprintf('  R3 (d) TI, Jw=1.44\n');
R4 = runcase(mdl, mws, withp(base,'quant_on',1));               fprintf('  R4 (e) conv + encoder\n');
R5 = runcase(mdl, mws, withp(base,'mode',2,'quant_on',1));      fprintf('  R5 (e) TI + encoder\n');

wa_list = [50 100 150 200 300 400 600 800];
B = cell(size(wa_list));
for i = 1:numel(wa_list)
    B{i} = runcase(mdl, mws, withp(base,'mode',2,'wa',wa_list(i),'act_on',1));
    fprintf('  bonus wa=%d\n', wa_list(i));
end
runset(mws, base);   % 기본값 복구

%% 2. (b)(c) 지표 ------------------------------------------------------------
fprintf('\n===== (b)(c) conventional vs transient-improved =====\n');
fprintf('%-10s | max|s| 0.5-1.5 | max|s| 1.5-3 | 오차@0.7[N·m] | 시간지연[ms] | 램프후 2%%정착[s] | 급락후 2%%정착[s] | Tm피크\n','케이스');
C = {R1, R2};  cn = {'conv', 'TI wa=100'};
for i = 1:2
    o = C{i};
    [te, e] = esterr(o);
    e07 = interp1(te, e, 0.7);
    fprintf('%-10s |     %6.3f     |    %6.3f    |    %7.1f    |    %6.1f    |      %6.4f      |      %6.4f      | %6.1f\n', ...
        cn{i}, maxs(o,0.5,1.5), maxs(o,1.5,3), e07, -e07/1500*1e3, ...
        settle(te, e, 0.7, 1.5, 6), settle(te, e, 1.5, 3, 3.6), max(abs(sq(o.Tm.Data))));
end
fprintf('해석값: conv 램프 정상오차 -%.0f N·m, TI -%.0f N·m(10 ms), TI 급락 2%% 정착 ln(50)/wa = %.1f ms\n', ...
    1.2*1500/293.88*1.2*20, 1500/100, log(50)/100*1e3);
fprintf('telescoped 검증 (R2): max|Td_hat - Td_hat_tele| = %.3e\n', ...
    max(abs(sq(R2.Td_hat.Data) - sq(R2.Td_hat_tele.Data))));

%% 3. (d) 지표 ---------------------------------------------------------------
fprintf('\n===== (d) Jw = 1.44, Jn = 1.2 (TI) =====\n');
[te3, e3] = esterr(R3);
W = [0.3 0.5; 0.9 1.4; 2.3 3.0];
for i = 1:size(W,1)
    m = mean(e3(te3 >= W(i,1) & te3 < W(i,2)));
    if W(i,1) < 2, an = 0.24*20; else, an = 0; end
    fprintf('구간 %.1f-%.1f s : 평균(Td_hat - Td) = %6.3f N·m  (해석 %.1f)\n', W(i,1), W(i,2), m, an);
end
fprintf('구동력 편향 (램프): %.1f N\n', 0.24*20/0.30);

%% 4. (e) 지표 ---------------------------------------------------------------
fprintf('\n===== (e) 4096 카운트 엔코더 차분 =====\n');
q = 2*pi/4096/1e-3;
fprintf('속도 양자화 단위 q = %.3f rad/s, q/sqrt(12) = %.3f, q/sqrt(6) = %.3f\n', q, q/sqrt(12), q/sqrt(6));
fprintf('%-10s | 이득[N·m/(rad/s)] | RMS(s 차이) | RMS(Tm 차이) | RMS(Td_hat 차이) | RMS(dTm)\n','케이스');
E = {R4, R5};  Cl = {R1, R2};  G = [1.2*20, 1.2*(20+100)];
for i = 1:2
    a = E{i};  c = Cl{i};
    idx = a.s.Time >= 0.2;
    ds  = sq(a.s.Data) - sq(c.s.Data);
    dT  = sq(a.Tm.Data) - sq(c.Tm.Data);
    dTh = sq(a.Td_hat.Data) - sq(c.Td_hat.Data);
    fprintf('%-10s |      %6.1f      |   %6.3f    |   %7.2f    |     %7.2f      | %7.2f\n', cn{i}, G(i), ...
        rmsv(ds(idx)), rmsv(dT(idx)), rmsv(dTh(idx)), rmsv(diff(sq(a.Tm.Data(idx)))));
end

%% 5. 보너스 -----------------------------------------------------------------
fprintf('\n===== 보너스: tau = 5 ms, wa sweep =====\n');
fprintf('%-6s | 급락후 언더슈트[%%] | 급락후 2%%정착[s] | std(s), t>2.2 | Tm피크[N·m]\n','wa');
[teR, eR] = esterr(R2);
fprintf('%-6s |      %6.1f        |     %6.4f       |   %8.5f    |  %6.1f   (지연 없음 기준)\n', 'ref', ...
    undershoot(R2), settle(teR, eR, 1.5, 3, 3.6), std(sq(R2.s.Data(R2.s.Time > 2.2))), max(abs(sq(R2.Tm.Data))));
for i = 1:numel(wa_list)
    o = B{i};  [tb, eb] = esterr(o);
    fprintf('%-6d |      %6.1f        |     %6.4f       |   %8.5f    |  %6.1f\n', wa_list(i), ...
        undershoot(o), settle(tb, eb, 1.5, 3, 3.6), std(sq(o.s.Data(o.s.Time > 2.2))), max(abs(sq(o.Tm.Data))));
end
diary off;

%% 6. 그림 -------------------------------------------------------------------
cC = [0 0.447 0.741];  cT = [0.466 0.674 0.188];

% (b) conventional 4종
f1 = figure('Name','P2 (b) conventional','Position',[50 50 900 900]);
t = R1.w.Time;
subplot(4,1,1); plot(t, sq(R1.w.Data), 'Color', cC); hold on; plot(t, wd_true(t), 'k--'); grid on
ylabel('\omega [rad/s]'); legend('\omega','\omega_d','Location','southeast')
subplot(4,1,2); plot(R1.s.Time, sq(R1.s.Data), 'Color', cC); grid on; ylabel('s [rad/s]')
subplot(4,1,3); plot(R1.Td.Time, sq(R1.Td.Data), 'k'); hold on
plot(R1.Td_hat.Time, sq(R1.Td_hat.Data), 'Color', cC); grid on
ylabel('[N\cdotm]'); legend('T_d','T_d hat','Location','northeast')
subplot(4,1,4); plot(R1.Tm.Time, sq(R1.Tm.Data), 'Color', cC); grid on
ylabel('T_m [N\cdotm]'); xlabel('t [s]')
savepng(f1, 'fig_P2b.png');

% (c) T_d hat 비교 + 확대
f2 = figure('Name','P2 (c) conv vs TI','Position',[80 80 900 750]);
subplot(2,2,[1 2]); hold on; grid on
plot(R1.Td.Time, sq(R1.Td.Data), 'k')
plot(R1.Td_hat.Time, sq(R1.Td_hat.Data), 'Color', cC)
plot(R2.Td_hat.Time, sq(R2.Td_hat.Data), 'Color', cT, 'LineWidth', 1.2)
ylabel('[N\cdotm]'); xlabel('t [s]'); legend('T_d','T_d hat conv','T_d hat TI \omega_a=100','Location','northeast')
subplot(2,2,3); hold on; grid on
plot(R1.Td.Time, sq(R1.Td.Data), 'k')
plot(R1.Td_hat.Time, sq(R1.Td_hat.Data), 'Color', cC)
plot(R2.Td_hat.Time, sq(R2.Td_hat.Data), 'Color', cT, 'LineWidth', 1.2)
xlim([0.45 1.1]); ylabel('[N\cdotm]'); xlabel('t [s]'); title('램프 구간')
subplot(2,2,4); hold on; grid on
plot(R1.Td.Time, sq(R1.Td.Data), 'k')
plot(R1.Td_hat.Time, sq(R1.Td_hat.Data), 'Color', cC)
plot(R2.Td_hat.Time, sq(R2.Td_hat.Data), 'Color', cT, 'LineWidth', 1.2)
xlim([1.45 2.1]); xlabel('t [s]'); title('마찰 급락')
savepng(f2, 'fig_P2c.png');

% (c) s 비교 (보고서에 넣을지 선택)
f3 = figure('Name','P2 (c) s 비교','Position',[110 110 900 300]);
plot(R1.s.Time, sq(R1.s.Data), 'Color', cC); hold on
plot(R2.s.Time, sq(R2.s.Data), 'Color', cT, 'LineWidth', 1.2); grid on
ylabel('s [rad/s]'); xlabel('t [s]'); legend('conv','TI \omega_a=100')
savepng(f3, 'fig_P2c_s.png');

% (d) 추정 편향
f4 = figure('Name','P2 (d) inertia','Position',[140 140 900 350]);
plot(te3, e3, 'Color', cT); hold on; grid on
plot([0 2], [4.8 4.8], 'r--', [2 3], [0 0], 'r--')
ylim([-20 20]); ylabel('T_d hat - T_d [N\cdotm]'); xlabel('t [s]')
legend('시뮬레이션 (J_w=1.44, J_{w,n}=1.2)','해석값 (J_w-J_{w,n})\omegȧ_d','Location','northeast')
savepng(f4, 'fig_P2d.png');

save('P2_results.mat','R1','R2','R3','R4','R5','B','wa_list');
fprintf('\n완료: P2_results.txt, P2_results.mat, fig_P2b/c/c_s/d.png 저장\n');


%% ======================= 함수 =======================
function o = runcase(mdl, mws, p)
    runset(mws, p);  o = sim(mdl);
end
function runset(mws, p)
    Jn = 1.2; lam = 20; Ts = 1e-3;
    assignin(mws,'mode',p.mode);     assignin(mws,'zeta',p.zeta);
    assignin(mws,'ka_conv', Jn^2*lam^2/(4*p.zeta^2));
    assignin(mws,'wa',p.wa);         assignin(mws,'ka',p.wa*Ts);
    assignin(mws,'Jw',p.Jw);         assignin(mws,'Jn',Jn);
    assignin(mws,'quant_on',p.quant_on);
    assignin(mws,'act_on',p.act_on); assignin(mws,'tau_act',5e-3);
end
function p = withp(p, varargin)
    for i = 1:2:numel(varargin), p.(varargin{i}) = varargin{i+1}; end
end
function v = sq(x), v = squeeze(x); end
function r = rmsv(x), r = sqrt(mean(x.^2)); end
function y = Td_true(t)
    y = min(max(1500*(t-0.5), 0), 300) - 180*(t >= 1.5);
end
function y = wd_true(t)
    y = min(20*t, 40);
end
function [t, e] = esterr(o)
    % 추정 오차 Td_hat - Td (컨트롤러 샘플 시각, Td는 해석식)
    t = o.Td_hat.Time;  e = sq(o.Td_hat.Data) - Td_true(t);
end
function m = maxs(o, t1, t2)
    t = o.s.Time; v = sq(o.s.Data);
    m = max(abs(v(t >= t1 & t < t2)));
end
function ts = settle(t, e, t0, t1, band)
    idx  = find(t >= t0 & t < t1);
    last = find(abs(e(idx)) > band, 1, 'last');
    if isempty(last), ts = 0; else, ts = t(idx(min(last+1, numel(idx)))) - t0; end
end
function u = undershoot(o)
    % 급락(300 -> 120) 후 Td_hat이 120 아래로 내려간 정도, 스텝 크기 180 대비 %
    t = o.Td_hat.Time; v = sq(o.Td_hat.Data);
    u = max(0, 120 - min(v(t >= 1.5 & t < 2.2))) / 180 * 100;
end
function savepng(f, name)
    try, exportgraphics(f, name, 'Resolution', 200); catch, saveas(f, name); end
end