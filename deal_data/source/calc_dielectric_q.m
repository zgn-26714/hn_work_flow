 clear
close all
clc

sizeoffont=16;
widthofline=2;
configPlot('FontSize',      sizeoffont,                ...
    'FontName',      'Helvetica', ...
    'Position',      [2 2   4.8*2*1.3 3.6*2*1.5],   ...
    'Units',        'centimeters',      ...
    'LineWidth',     2,                 ...
    'AxesLineWidth', 1.2);
%%
all_mol = {'water','pc','dmso','acn','mt','ac','dme'};
for my_idx=1:length(all_mol)
    %% calc J
    mol=all_mol{my_idx};
    
    disp(mol);
    all_data = load(sprintf('./fft_data/%s_dphiANDdieletric.mat',mol));
    fprintf('  数据加载完成，正在计算J_omega...\n');
%     f = all_data.dielectric.f;

    f_max = 1e17;  % Hz
    n_points = 10000;
    f_min = 1e5;
%     clear f
    f = logspace(log10(f_min), log10(f_max), n_points);
    f = [0, f];


    omega = 2*pi*f;
    omega=omega(:);
    %calc eps
    A_list = all_data.dielectric.A_list;
    tau_list = all_data.dielectric.tau_list;
    eps_static=all_data.dielectric.eps_static;
    epsilon_inf=all_data.dielectric.epsilon_inf;
    num   = (A_list(:).*tau_list(:)).';     % 1×Nexp
    denom = 1 + 1i*(omega * tau_list(:).'); % nω×Nexp
    S = sum(num ./ denom, 2);               % nω×1
    delta_eps = eps_static - epsilon_inf;
    S=S(:);
    eps_an    = epsilon_inf + delta_eps * (1 - 1i*omega.*S);


    %calc dphi

    fprintf('\n加载onfly数据...\n');
    onfly_data = load(sprintf('C:\\Users\\Administrator\\Desktop\\mycode\\JEModel\\Qresult\\notsmooth\\notfit\\%s298k2v0psIE.mat',mol));
    t_onfly = linspace(onfly_data.onfly_dt, size(onfly_data.JEQ,2)*onfly_data.onfly_dt, size(onfly_data.JEQ,2));

    E_a = all_data.fit_params.a;
    tau_s = all_data.fit_params.tau_s;

    E_omega = zeros(size(omega));
    for i=1:length(E_a)
        E_omega = E_omega + E_a(i)*tau_s(i)./(1 + 1i*omega*tau_s(i));
    end
%     E_omega = E_omega + (0.2 - sum(E_a))./(1i * omega);
    
    
%     eps_an=eps_an(:);
    E_omega=E_omega(:);
    J_omega= omega .* (1i) .* (eps_an-1) .* E_omega * 8.854187817e-12;
    J_omega = J_omega * 6.2415e-3; %% unit e / ps /nm^2
    [J_omega, omega]=restore_neg(J_omega, omega); %% 拓展全频域
    fprintf('  全频域扩展完成\n');
    fprintf('\n进行傅里叶逆变换...\n');
    [t_current, J_t]=inverse_fourier_numeric(J_omega, omega, onfly_data.onfly_dt, size(onfly_data.JEQ,2)*onfly_data.onfly_dt);  %% 傅里叶逆变换
    figure
    plot(t_current*1e12, -J_t,'b');hold on
    plot(t_onfly, mean(onfly_data.J_div{1}(300:700,:)),'r--')
    xlabel('{\itt} (ps)')
    ylabel('{\itq} (e/ps/nm^2)')
    legend('dieletric', 'JEmodel');
    set(gca,'Units','centimeters','Position',[4 2.5 7 6])
    xlim([0 size(onfly_data.JEQ,2)*onfly_data.onfly_dt])
    title('current')

    figure
    plot(omega,abs(J_omega));hold on
    %% calc q
    fprintf('\n计算q_omega...\n');
    qf_max = 1e16;  % Hz
    n_points = 10000;
    qf_min = 1e6;
    qf = logspace(log10(qf_min), log10(qf_max), n_points);
    qf = [0.001, qf];
    
    omegaOFq = 2 * pi * qf;
    fprintf('  生成频率点: %d 个点，从 %.2e Hz 到 %.2e Hz\n', n_points, qf_min, qf_max);

    q_omega = get_q_omega(J_omega, all_data.fit_params, omega, omegaOFq);%% unit eV/ps/nm^3

    fprintf('  q_omega计算完成，正在扩展到全频域...\n');
    [q_omega, omegaOFq]=restore_neg(q_omega, omegaOFq);  %% 拓展全频域

    fprintf('\n进行傅里叶逆变换...\n');
    [t_heat, q_t]=inverse_fourier_numeric(q_omega, omegaOFq, onfly_data.onfly_dt, size(onfly_data.JEQ,2)*onfly_data.onfly_dt);  %% 傅里叶逆变换

    figure;
    plot(t_heat*1e12,q_t)
    xlim([0 size(onfly_data.JEQ,2)*onfly_data.onfly_dt])
    xlabel('{\itt} (ps)')
    ylabel('{\itq} (eV/ps/nm^3)')
    hold on
    plot(t_onfly, mean(onfly_data.JEQ_div{1}(300:700,:))/96.488,'r--')
    legend('dieletric', 'JEmodel');
    set(gca,'Units','centimeters','Position',[4 2.5 7 6])




    fprintf('\n=============== 所有计算完成 ===============\n');

    savefile=sprintf('./Qresult/%s_qomega.mat',mol);
    save(savefile, 'omega','J_omega', 'q_omega','omegaOFq','J_t','q_t','t_current','t_heat','E_omega','eps_an');

    figure;
    plot(omegaOFq/2/pi, imag(q_omega))
    xlabel('{\itf} (Hz)')
    ylabel('{\itq} (\omega)')
    set(gca,'Units','centimeters','Position',[4 2.5 7 6])
    set(gca, 'XScale', 'log');
    xlim([qf_min qf_max])
    xticks([qf_min  qf_max])
end

function q_omega = get_q_omega(J_omega, fit_para, omega,omegaOFq)
    omega = omega(:);
    J_omega = J_omega(:);
    N = length(omega);
    E_a = fit_para.a;
    tau_s = fit_para.tau_s;
    
    % 积分权重（梯形法则）
    deltaOmega = zeros(N,1);
    deltaOmega(1) = (omega(2)-omega(1))/2;
    for m = 2:N-1
        deltaOmega(m) = (omega(m+1)-omega(m-1))/2;
    end
    deltaOmega(N) = (omega(N)-omega(N-1))/2;
    
    q_omega = zeros(length(omegaOFq),1);
    
    h = waitbar(0, '正在进行逆傅里叶变换...');
    hn_N = length(omegaOFq);   % 循环次数
    
    
    for k = 1:length(omegaOFq)
        waitbar(k/hn_N, h, sprintf('进度：%.1f %%', k/hn_N*100));

        % 计算卷积：Q(ω_k) = ∫ J(ω') E(ω_k - ω') dω'
        omega_diff = omegaOFq(k) - omega;  % ω_k - ω'
        
        % 用解析式计算 E(ω_k - ω')
%         omega_diff(omega_diff == 0) = nan;
        E_diff = zeros(size(omega));
        for i=1:length(E_a)
            E_diff = E_diff + E_a(i)*tau_s(i)./(1 + 1i*omega_diff*tau_s(i));
        end
%         E_diff = E_diff +  (0.2 - sum(E_a))./(1i * omega_diff);
%         E_diff(isnan(E_diff)) = 0;
%         J_omega(isnan(J_omega)) = 0;
%         tmp = J_omega .* E_diff .* deltaOmega;
%         tmp(isnan(tmp) | isinf(tmp)) = 0;
%         idx = (length(tmp)-1)/2;
%         front_tmp = sum(tmp(1:idx));
%         behind_tmp = sum(tmp(idx+2:end));
        % 数值积分
        q_omega(k) = sum(J_omega .* E_diff .* deltaOmega) / (2*pi);
    end
    close(h);  
end
                
function [t, x] = inverse_fourier_numeric(X, omega, onfly_dt, all_t, method)
% 数值积分法逆傅里叶变换
% 输入：
%   X - 频域信号（复数）
%   omega - 角频率向量（rad/s）
%   method - 积分方法：'trapz'（默认）, 'simpson', 'quad'
% 输出：
%   t - 时间向量
%   x - 时域信号

if nargin < 5
    method = 'trapz';
end

% 确定时间向量
t = linspace(onfly_dt, all_t, 10000);
% 预分配时域信号
t = t *1e-12;
x = zeros(size(t));

h = waitbar(0, '正在进行逆傅里叶变换...');
N = length(t);   % 循环次数

% 根据选择的积分方法计算
switch lower(method)
    case 'trapz'
        % 梯形法数值积分
        for k = 1:length(t)
            integrand = X .* exp(1i*omega*t(k));
            x(k) = trapz(omega, integrand) / (2*pi);
             waitbar(k/N, h, sprintf('进度：%.1f %%', k/N*100));
        end
        
    case 'simpson'
        % Simpson法数值积分（更精确）
        for k = 1:length(t)
            integrand = X .* exp(1i*omega*t(k));
            x(k) = simpson_integration(omega, integrand) / (2*pi);
             waitbar(k/N, h, sprintf('进度：%.1f %%', k/N*100));
        end
        
    case 'quad'
        % 自适应积分（精度最高，但最慢）
        warning('自适应积分较慢，请耐心等待...');
        for k = 1:length(t)
            x(k) = integral(@(w) integrand_func(w, X, omega, t(k)), ...
                           min(omega), max(omega)) / (2*pi);
                        waitbar(k/N, h, sprintf('进度：%.1f %%', k/N*100));
        end
        
    otherwise
        error('不支持的积分方法');
end
close(h);  
% 取实部（理论上结果应为实数）
x = real(x);

end

function I = simpson_integration(x, y)
% Simpson 1/3法则数值积分
n = length(x);
if mod(n,2) == 0
    % 偶数个点，使用复合Simpson
    h = (x(end) - x(1)) / (n-1);
    I = (h/3) * (y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end));
else
    % 奇数个点，直接使用Simpson
    I = (x(end)-x(1))/6 * (y(1) + 4*y((1+end)/2) + y(end));
end
end

function y = integrand_func(w, X, omega, t)
% 插值计算被积函数
X_interp = interp1(omega, X, w, 'spline', 0);
y = X_interp .* exp(1i*w*t);
end

function [J_omega_full, omega_full] = restore_neg(J_omega, omega)
% RESTORE_NEG - 为正频域数据补全负频率部分，使信号满足厄米对称
%
% 输入：
%   omega   : 原始频率数组（通常为正频率, 单位可以是 rad/s, Hz等）
%   J_omega : 对应的频域函数值（可以是复数）
%
% 输出：
%   omega_full : 包含负频率 + 正频率的完整频率轴（从负到正）
%   J_omega_full : 完整频谱，满足 J(-w)=conj(J(w))
%
% 注意：
%   1. 频率轴 omega 必须从小到大（如果不是，会自动排序）
%   2. 应确保 omega(1) = 0（如果没有 0 点也没关系）
%

    % --- 确保频率从小到大排序 ---
    [omega, idx] = sort(omega(:));
    J_omega = J_omega(idx);

    % --- 判断是否包含 0 频率 ---
    hasZero = any(abs(omega) < 1e-15);

    % --- 构建负频率 ---
    if hasZero
        % 如果存在 0 Hz：负频率部分为所有正频率的镜像（不包含0）
        omega_neg = -flipud(omega(2:end));
        J_neg = conj(flipud(J_omega(2:end)));  
        
        omega_full = [omega_neg; omega];
        J_omega_full = [J_neg; J_omega];
    else
        % 如果没有 0 点：全部点都要镜像
        omega_neg = -flipud(omega);
        J_neg = conj(flipud(J_omega));

        omega_full = [omega_neg; omega];
        J_omega_full = [J_neg; J_omega];
    end
end



function [t, Jt] = analytic_Jt(A_list, tau_list, Ea, tau_s, eps_static, epsilon_inf, t_max, n_t)
% analytic_Jt  计算解析的 J(t)
% 输入:
%   A_list, tau_list        : 电介质分项（列向量 or 行向量）
%   Ea, tau_s               : 电场分项（列向量 or 行向量），Ea 对应 tau_s
%   eps_static, epsilon_inf : 标量
%   t_max                   : 最大时间（单位同 tau 的单位）；例如 200 （ps）
%   n_t                     : 时间点数，例如 1000
%
% 输出:
%   t  : 时间向量（同 tau 单位）
%   Jt : J(t) （单位与脚本中换算一致： e/ps/nm^2）

% constants
eps0 = 8.854187817e-12;      % F/m
conv = 6.2415e-3;            % 乘子（你原代码）
C = eps0 * conv;

% ensure column vectors
A_list = A_list(:);
tau_list = tau_list(:);
Ea = Ea(:);
tau_s = tau_s(:);

Delta_eps = eps_static - epsilon_inf;

% time vector：假设 tau 单位为 ps（与原代码一致），生成 t 从 0 到 t_max（ps）
t = linspace(0, t_max, n_t).';   % column vector, units: ps

% preallocate
Jt = zeros(size(t));

% first term J^(A): (eps_static - 1) * d/dt E(t)
% E(t) = sum_k Ea_k * exp(-t/tau_s_k)
% d/dt E(t) = sum_k Ea_k * (-1/tau_s_k) * exp(-t/tau_s_k)
termA = zeros(size(t));
for kk = 1:length(Ea)
    termA = termA + Ea(kk) * (-1./tau_s(kk)) .* exp(-t./tau_s(kk));
end
termA = (eps_static - 1) * termA;   % note sign included

% second term: B = Delta_eps * omega^2 * (sum_j G_j) * E_omega  -> time domain: -Delta_eps * d^2/dt^2 h(t)
% where h(t) = sum_{j,k} A_j Ea_k * (tau_j tau_s_k)/(tau_s_k - tau_j) * ( e^{-t/tau_s_k} - e^{-t/tau_j} )
tol = 1e-12; % tolerance for equality of taus (adjust if tau in other units)
termB = zeros(size(t));

for j = 1:length(A_list)
    Aj = A_list(j);
    tj = tau_list(j);
    for k = 1:length(Ea)
        Ek = Ea(k);
        ts = tau_s(k);
        pref = Aj * Ek;
        if abs(ts - tj) > 1e-9 % distinct taus
            coeff = (tj * ts) / (ts - tj);
            % h_jk(t) = pref * coeff * (e^{-t/ts} - e^{-t/tj})
            % compute d^2/dt^2 h_jk = pref*coeff*( 1/ts^2 e^{-t/ts} - 1/tj^2 e^{-t/tj} )
            d2h = pref * coeff .* ( (1/(ts^2)) .* exp(-t./ts) - (1/(tj^2)) .* exp(-t./tj) );
        else
            % ts == tj (within tol): use limit -> h_jk(t) -> pref * t * exp(-t/tau)
            tau = (ts + tj)/2;
            % h = pref * t * exp(-t/tau)
            % d2h/dt2 = pref * d2( t e^{-t/tau} )/dt2
            % first derivative: d/dt = pref * ( e^{-t/tau} - (t/tau) e^{-t/tau} )
            % second derivative: d2/dt2 = pref * ( - (1/tau) e^{-t/tau} - (1/tau) e^{-t/tau} + (t/tau^2) e^{-t/tau} )
            % simplifying: pref * e^{-t/tau} * ( t/tau^2 - 2/tau )
            d2h = pref .* exp(-t./tau) .* ( (t./(tau^2)) - 2./tau );
        end
        % termB accumulates -Delta_eps * d2h
        termB = termB - Delta_eps * d2h;
    end
end

% total (before global constant)
y = termA + termB;

% multiply by conversion constant C
Jt = C .* y;

% enforce causality (t<0 zero) — t vector already >=0
end


function q_omega_an = cal_q_an(tau_1, tau_2, C, omega)
% C=delta_eps*E_a^2*tau_s^2/(tau_s-tau_list);

C_q = C*tau_2^2/(tau_2-tau_1);
left = 1./(tau_1 + tau_2 + 1i.*omega.*tau_1.*tau_2);
right = 1./ (2*tau_2+1i .* omega .* tau_2^2);
q_omega_an = C_q.*(left-right);
q_omega_an = q_omega_an *8.854187817e-12 *6.2415e-3;
end
