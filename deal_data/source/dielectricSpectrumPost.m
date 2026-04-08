clear
% close all
clc

sizeoffont=16;
widthofline=2;
configPlot('FontSize',      sizeoffont,                ...
    'FontName',      'Helvetica', ...
    'Position',      [2 2   4.8*2*1.3 3.6*2*1.5],   ...
    'Units',        'centimeters',      ...
    'LineWidth',     2,                 ...
    'AxesLineWidth', 1.2);
%% ------------------ 用户区：参数 ------------------
all_mol = {'water','pc','dmso','acn','mt','ac','dme'};

V_all=[88.0852 88.1943 86.1342 87.1448 87.5229 86.845 84.0957];
for my_idx=1:length(all_mol)
    clearvars -except all_mol V_all my_idx
    mol=all_mol{my_idx};
    files.Mtot   =  sprintf('./data/%s_Mtot.xvg',mol);   % 必需：总偶极矩时间序列：t(ps) Mx(My,Mz) (Debye)
    files.epsw   =  sprintf('./data/%s_epsw.xvg',mol);   % 可选：GROMACS 频域 ε' ε'' 对比
    files.dipcorr=  sprintf('./data/%s_dipcorr.xvg',mol);% 可选：GROMACS 归一化相关函数
    files.epsilon=  sprintf('./data/%s_epsilon.xvg',mol);% 可选：GROMACS ε(t)
    files.deriv  =  sprintf('./data/%s_deriv.xvg',mol);  % 可选：GROMACS dε/dt

    % 频率轴
    f_max = 1e30;  % Hz
    n_points = 500000;
    f_min = 1e3;
    f = logspace(log10(f_min), log10(f_max), n_points);
    f = [0, f];

    T       = 298;                 % 温度 (K)
    L_nm    = 5.319;               % 立方盒边长 (nm) -> 体积 L^3
    cutN    = 7000;                 % 可选：只用前 N 行（避免巨量数据演示慢）

    DO_FIT  = true;                % 是否做多指数拟合
    N_EXP   = 2;                   % 多指数个数 (1,2,3…)
    SPECTRUM_METHOD = 'both';      % 'fft' | 'analytic' | 'both'

    epsilon_inf = 1.0;             % 高频极限 ε∞，如需可改

    % 拟合选项
    fit_options = optimset('MaxIter',5000,'MaxFunEvals',5000,'TolX',1e-7,'TolFun',1e-7,'Display','off');
    plot_ps_xlim = 500;            % 时域图 x 轴上限 (ps)

    %% ------------------ 常量与单位 ------------------
    eps0        = 8.854187817e-12; % F/m
    kB          = 1.380649e-23;    % J/K
    D_to_Cm     = 3.33564e-30;     % 1 Debye = 3.33564e-30 C·m
    nm3_to_m3   = 1e-27;           % 1 nm^3 = 1e-27 m^3

    V_nm3 = V_all(my_idx);



    V = V_nm3 * nm3_to_m3;

    %% ------------------ 数据读取 ------------------
    Mdata = xvgread_simple(files.Mtot);
    if ~isempty(cutN) && size(Mdata,1) > cutN, Mdata = Mdata(1:cutN,:); end

    epsw    = safe_read_xvg(files.epsw);
    dipcorr = safe_read_xvg(files.dipcorr);
    epsilon = safe_read_xvg(files.epsilon);
    deriv   = safe_read_xvg(files.deriv);

    time_ps = Mdata(:,1);
    Mx = Mdata(:,2) * D_to_Cm;
    My = Mdata(:,3) * D_to_Cm;
    Mz = Mdata(:,4) * D_to_Cm;

    dt_ps = time_ps(2) - time_ps(1);
    dt    = dt_ps * 1e-12;
    N     = length(time_ps);
    t_corr = (0:N-1)' * dt;   % s

    % 偶极统计
    M_total = [Mx, My, Mz];
    M_mean  = mean(M_total,1);
    M_fluct = M_total - M_mean;
    M_sq = sum(M_total.^2,2);
    M2_mean = mean(M_sq);
    M_mean_sq = sum(M_mean.^2);
    fluct = M2_mean - M_mean_sq;

    % 预因子与涨落 εs
    prefactor  = 1/(3*eps0*V*kB*T);
    eps_static = 1 + prefactor*(M2_mean - M_mean_sq);

    %% ------------------ 自相关函数（标量形式） ------------------
    % 定义：C(t) = <ΔM(0)·ΔM(t)>，并归一化 C(0)=1
    C = autocorr_dot(M_fluct); %自相关函数的计算
    C0 = C(1);
%     C_norm = C / C0;
    C_norm = C;

    %% ------------------ 拟合：多指数（可选） ------------------
    A_list = 1; tau_list = inf; C_fit_norm = C_norm; % 默认占位
    if DO_FIT
        [A_list, tau_list, C_fit_norm] = fit_multi_exp(C_norm, t_corr, N_EXP, dt, fit_options);
    end

    %% ------------------ 频域：FFT/解析（可选） ------------------
    omega = 2*pi*f;
    omega=omega(:);
    % 结果容器
    spec = struct();

    Nfft = 2^nextpow2(4*N);
    freq_full = (0:(Nfft-1))' / (Nfft*dt);
    idx_pos = 1:floor(Nfft/2);
    f_fft = freq_full(idx_pos);
    omega_fft = 2*pi.*f_fft;
    
    % 1) FFT 路径
    if any(strcmpi(SPECTRUM_METHOD, {'fft','both'}))
        C_for_fft = C_fit_norm.';               % 若拟合，则对拟合曲线做 FFT 更平滑；你也可以换 C_norm
        Cw_full   = fft(C_for_fft, Nfft) * dt;
        Cw        = Cw_full(idx_pos);
        delta_eps = eps_static - epsilon_inf;
        eps_fft   = epsilon_inf + delta_eps * (1 - 1i*omega_fft.*Cw);
        spec.fft.real = real(eps_fft);
        spec.fft.imag = -imag(eps_fft);         % 使损耗为正
    end

    % 2) 解析 路径（多指数解析）
    if DO_FIT && any(strcmpi(SPECTRUM_METHOD, {'analytic','both'}))
        % S(ω) = Σ_i A_i τ_i / (1 + i ω τ_i)
        num   = (A_list(:).*tau_list(:)).';     % 1×Nexp
        denom = 1 + 1i*(omega * tau_list(:).'); % nω×Nexp
        S = sum(num ./ denom, 2);               % nω×1
        delta_eps = eps_static - epsilon_inf;
        eps_an    = epsilon_inf + delta_eps * (1 - 1i*omega.*S);
        spec.an.real = real(eps_an);
        spec.an.imag = -imag(eps_an);
    end

    %% ------------------ 特征量（各自方法独立给出） ------------------
    stat = struct();
    if isfield(spec,'fft')
        [stat.fft.im_max, imax] = max(spec.fft.imag);
        stat.fft.f_max  = f(imax);
        stat.fft.tau_ps = 1/(2*pi*stat.fft.f_max)*1e12;
    end
    if isfield(spec,'an')
        [stat.an.im_max, imax] = max(spec.an.imag);
        stat.an.f_max  = f(imax);
        stat.an.tau_ps = 1/(2*pi*stat.an.f_max)*1e12;
    end

    %% ------------------ 绘图（精简版） ------------------
    % 1. 时域相关 (线型 & 拟合)
    figure('Name','C(t)');plot(t_corr*1e12, C_norm, 'b-', 'LineWidth', 1.8); hold on;
    if DO_FIT, plot(t_corr*1e12, C_fit_norm, 'r--', 'LineWidth', 1.6); end
    xlabel('t (ps)'); ylabel('C(t) (norm)'); title('Autocorr');
    grid on; 
    xlim([0 min(plot_ps_xlim, t_corr(end)*1e12)]);
    legend({'raw','fit'},'Location','best');

    figure('Name','\epsilon'' ');
    if isfield(spec,'fft'), semilogx(f_fft, spec.fft.real, 'b-', 'LineWidth', 1.8); hold on; end
    if isfield(spec,'an'),  semilogx(f, spec.an.real,  'k--','LineWidth', 1.8); end
    xlabel('f (Hz)'); ylabel("\epsilon'"); title('实部');
    % grid on; 
    legend(find_legend(spec,'real'),'Location','best','Interpreter','tex');

    figure('Name','\epsilon'''' ');
    if isfield(spec,'fft'), semilogx(f_fft, spec.fft.imag, 'r-', 'LineWidth', 1.8); hold on; end
    if isfield(spec,'an'),  semilogx(f, spec.an.imag,  'k--','LineWidth', 1.8); end
    xlabel('f (Hz)'); ylabel("\epsilon''"); title('虚部');
    % grid on; 
    legend(find_legend(spec,'imag'),'Location','best');

    savefile=sprintf('savefile/%s_epsilon.mat',mol);
    real=spec.an.real;
    imag=spec.an.imag;
    save(savefile,'f','real','imag','t_corr','C_norm','C_fit_norm','spec','f_fft','epsw');



    % 4. Cole-Cole
    figure('Name','Cole-Cole');
    if isfield(spec,'fft'), plot(spec.fft.real, spec.fft.imag, 'b-','LineWidth',1.8); hold on; end
    if isfield(spec,'an'),  plot(spec.an.real,  spec.an.imag,  'k--','LineWidth',1.8); end
    xlabel("\epsilon'"); ylabel("\epsilon''"); title('Cole-Cole'); axis equal; 
    % ...grid on;
    legend(find_legend(spec,'cole'),'Location','best');

    figure('Name','Compare GMX');
    has_g = ~isempty(epsw);
    if has_g
        semilogx(epsw(:,1)*1e9, epsw(:,2), 'b:', 'LineWidth',1.2); hold on;
        semilogx(epsw(:,1)*1e9, -epsw(:,3), 'r:', 'LineWidth',1.2);
    end
    if isfield(spec,'fft')
        semilogx(f_fft, spec.fft.real, 'b-', 'LineWidth',1.6); hold on;
        semilogx(f_fft, spec.fft.imag, 'r-', 'LineWidth',1.6);
    end
    if isfield(spec,'an')
        semilogx(f, spec.an.real, 'k--', 'LineWidth',1.6);
        semilogx(f, spec.an.imag, 'k-.', 'LineWidth',1.6);
    end
    xlabel('f (Hz)'); ylabel('\epsilon(\omega)'); title('Compare (if GROMACS present)'); 
    % grid on;
    lgd = {};
    if has_g, lgd = [lgd, "\epsilon' GROMACS","\epsilon'' GROMACS"]; end
    if isfield(spec,'fft'), lgd = [lgd, "\epsilon' FFT","\epsilon'' FFT"]; end
    if isfield(spec,'an'),  lgd = [lgd, "\epsilon' Analytic","\epsilon'' Analytic"]; end
    legend(lgd,'Location','best');

    %% ------------------ 输出关键信息到命令行 ------------------
    fprintf('\n=== Dielectric Spectrum Analysis Results ===\n');
    fprintf('Temperature: T = %.1f K\n', T);
    fprintf('Simulation box volume: V = %.3f nm^3\n', V_nm3);
    fprintf('Static dielectric constant (fluctuation): eps_s = %.3f\n', eps_static);
    fprintf('High-frequency dielectric constant: eps_inf = %.3f\n', epsilon_inf);
    fprintf('Multi-exponential fit terms: N = %d\n', N_EXP);

    if isfield(stat,'fft')
        fprintf('\n--- FFT Method Results ---\n');
        fprintf('Loss peak frequency: f_max = %.2f GHz\n', stat.fft.f_max/1e9);
        fprintf('Relaxation time: tau = %.2f ps\n', stat.fft.tau_ps);
        fprintf('Maximum loss value: eps''''_max = %.3f\n', stat.fft.im_max);
    end

    if isfield(stat,'an')
        fprintf('\n--- Analytic Method Results ---\n');
        fprintf('Loss peak frequency: f_max = %.2f GHz\n', stat.an.f_max/1e9);
        fprintf('Relaxation time: tau = %.2f ps\n', stat.an.tau_ps);
        fprintf('Maximum loss value: eps''''_max = %.3f\n', stat.an.im_max);
    end

    % 如果有GROMACS数据
    if ~isempty(epsw)
        % 找到GROMACS数据的峰值
        [gmx_im_max, gmx_idx] = max(-epsw(:,3));  % 注意负号，因为GROMACS输出的是负值
        gmx_f_max = epsw(gmx_idx,1) * 1e9; % 转换为GHz

        % 检查频率是否合理
        if gmx_f_max > 1e6  % 如果大于1THz，说明单位可能有问题
            gmx_f_max = epsw(gmx_idx,1);  % 直接用原始值
        end

        gmx_tau_ps = 1/(2*pi*gmx_f_max*1e9) * 1e12; % 转换为ps

        fprintf('\n--- GROMACS Reference Results ---\n');
        fprintf('Loss peak frequency: f_max = %.2f GHz\n', gmx_f_max);
        fprintf('Relaxation time: tau = %.2f ps\n', gmx_tau_ps);
        fprintf('Maximum loss value: eps''''_max = %.3f\n', gmx_im_max);
    end

    % 输出拟合参数
    if DO_FIT && N_EXP >= 1
        fprintf('\n--- Multi-exponential Fit Parameters ---\n');
        for i = 1:length(A_list)
            fprintf('Component %d: weight A_%d = %.3f, relaxation time tau_%d = %.2f ps\n', ...
                i, i, A_list(i), i, tau_list(i)*1e12);
        end
        fprintf('Fit quality check: C_fit(0) = %.6f\n', C_fit_norm(1));
    end

    fprintf('\n=== Analysis Complete ===\n');

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%------------------ 第二部分：电场的傅里叶变换 ------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 读取数据
    data = load(sprintf('C:\\Users\\Administrator\\Desktop\\mycode\\JEModel\\Electric_field\\%s_elecField2V0ps.mat', mol));
    t = data.t(:);  % 时间单位：ps
    E_field = mean(data.elec_filed(300:700,:)).';
    E_field = E_field(:);

    %% =============================
    %       可选指数数量 N
    % ==============================
    N_exp = 2;   % <<<<<<<<<< 你可以改成 1, 2, 3, 4 ... 任意
    % ==============================

    %% 自动构建多指数模型字符串

    E0 = 0.2; 
    % 构建新模型：E(t) = E0 - sum( a_k * (1 - exp(-x/tau_k)) )
    model_terms = strings(1, N_exp);
    param_names = cell(1, 2*N_exp);  % 存所有参数名

    for k = 1:N_exp
        model_terms(k) = sprintf('a%d*(1 - exp(-x/tau%d))', k, k);
        param_names{2*k-1} = sprintf('a%d', k);
        param_names{2*k}   = sprintf('tau%d', k);
    end

    model_str = sprintf('%f - (%s)', E0, strjoin(model_terms, ' + '));

    % 构建 fittype
    double_exp_model = fittype(model_str, ...
        'independent','x', 'dependent','y', ...
        'coefficients',param_names);


    initial_params = ones(1,2*N_exp);
    lower_bounds = [ zeros(1,N_exp), 0.01*ones(1,N_exp) ];
    upper_bounds = Inf(1,2*N_exp);

    fit_options = fitoptions('Method','NonlinearLeastSquares',...
        'StartPoint',initial_params,...
        'Lower',lower_bounds,...
        'Upper',upper_bounds);

    [fit_result, gof] = fit(t, E_field, double_exp_model, fit_options);


    disp('========== 拟合参数 ==========');
    for k = 1:N_exp
        fprintf('a%d   = %.6e\n', k, fit_result.(sprintf('a%d',k)));
        fprintf('tau%d = %.6f ps\n', k, fit_result.(sprintf('tau%d',k)));
    end

    disp(gof);

    fprintf('=====================================\n\n');

    %% 3. Plot time domain fitting results
    figure;
    plot(t, E_field, 'b', 'LineWidth', 2, 'DisplayName', 'Original Data');
    hold on;
    plot([0 t'], fit_result([0 t']), 'r-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
    xlabel('{\itt} (ps)');
    ylabel('E (V/nm)');
    legend('Location', 'best');
    title('Time Domain Fitting');

    %% 4. Fourier transform of analytical expression

    a = zeros(1, N_exp);
    tau_s = zeros(1, N_exp);

    % 拟合参数 → 秒
    for k = 1:N_exp
        a(k) = fit_result.(sprintf('a%d', k));
        tau_ps = fit_result.(sprintf('tau%d', k));
        tau_s(k) = tau_ps * 1e-12;
    end

    FT_analytical = zeros(size(f));

    % 核心公式： a_k * tau_k / (1 + i 2π f tau_k)
    for k = 1:N_exp
        FT_analytical = FT_analytical + ...
            a(k) * tau_s(k) ./ (1 + 1i * 2*pi*f*tau_s(k));
    end
%     FT_analytical = FT_analytical + 

    amplitude = abs(FT_analytical);
    phase = angle(FT_analytical);

    %% 5. Plot frequency domain results
    figure
    plot(f, amplitude, 'b-', 'LineWidth', 2);
    hold on
    ylabel('{\itE} (\omega)');
    xlabel('Frequency (Hz)');
    xlim([1e7 1e13])
    set(gca, 'XScale', 'log');
    title('Analytical Fourier Transform');


    yyaxis right
    xlim([1e5 1e15])
    xticks([1e7 1e10 1e13])
    plot(f,spec.an.imag,'r')
    ylabel('\epsilon (\omega)')
    title(sprintf('%s',mol))

    fit_params = struct('a', a, 'tau_s', tau_s, ...
        'rsquare', gof.rsquare, 'rmse', gof.rmse);

    frequency_domain = struct('frequency', f, 'amplitude', amplitude, ...
        'phase', phase, 'complex', FT_analytical);

    dielectric = struct( ...
        'f',            f, ...
        'eps_an',       eps_an,...
        'real',         spec.an.real, ...
        'imag',         spec.an.imag, ...
        't_corr',       t_corr, ...
        'C_norm',       C_norm, ...
        'C_fit_norm',   C_fit_norm, ...
        'A_list',       A_list,...
        'tau_list',     tau_list,...
        'eps_static',   eps_static,...
        'epsilon_inf',  epsilon_inf...
    );

    savefile=sprintf('./fft_data/%s_dphiANDdieletric.mat',mol);
    save(savefile, 'fit_params', 'frequency_domain', 'dielectric');
end




%% ====================== 内部函数 ======================
function C = autocorr_dot(M)
% 标量自相关：C(t) = <ΔM(0)·ΔM(t)>
Nloc = size(M,1);
C = zeros(Nloc,1);
for tau = 1:Nloc
    n_terms = Nloc - tau + 1;
    dotp = sum(M(1:n_terms,:).*M(tau:Nloc,:), 2);
    C(tau) = mean(dotp);
end
end

function [A_list, tau_list, Cfit] = fit_multi_exp(Cnorm, t, nExp, dt, opts)
% 拟合：C_fit(t) = Σ A_i exp(-t/τ_i), ΣA_i=1, A_i>=0, τ_i>0
y = Cnorm(:); x = t(:);
% 初值（对数均分取 τ 初猜）
Tmax = max(x);
tau0 = log( linspace(10*dt, max(Tmax/4, 20*dt), nExp) )';
a0   = zeros(nExp,1)*Cnorm(1);                  % 幅度logits 初值
p0   = [a0; tau0];                     % [a, log(τ)]
fun  = @(p) sum( (model_exp_mix(p, x)' - y').^2 );
p_opt = fminsearch(fun, p0, opts);
[A_list, tau_list] = decode_params(p_opt, nExp);
Cfit = model_exp_mix(p_opt, x).';
end

function Cmodel = model_exp_mix(p, t)
n = numel(p)/2; a = p(1:n); b = p(n+1:end);
ea = exp(a - max(a));
A = ea / sum(ea);
tau = exp(b);
Cmodel = zeros(size(t));
for k = 1:n
    Cmodel = Cmodel + A(k)*exp(-t/tau(k));
end
end

function [A, tau] = decode_params(p, n)
a = p(1:n); b = p(n+1:end);
ea = exp(a - max(a));
A  = ea / sum(ea);
tau= exp(b);
end

function data = safe_read_xvg(fname)
if exist(fname,'file'), data = xvgread_simple(fname); else, data = []; end
end

function M = xvgread_simple(fname)
% 简单 xvg 读取：跳过以 '@' 或 '#' 开头的行，返回数值矩阵
fid = fopen(fname,'r');
assert(fid>0, 'Cannot open file: %s', fname);
C = textscan(fid, '%s', 'Delimiter','\n','Whitespace','','ReturnOnError',false);
fclose(fid);
lines = C{1};
vals = [];
for i=1:numel(lines)
    s = strtrim(lines{i});
    if isempty(s) || startsWith(s,'@') || startsWith(s,'#'), continue; end
    nums = sscanf(s, '%f');
    if ~isempty(nums), vals = [vals; nums(:)']; end %#ok<AGROW>
end
M = vals;
end

function lg = find_legend(spec, mode)
lg = {};
switch mode
    case 'real'
        if isfield(spec,'fft'), lg{end+1} = '\epsilon'' FFT'; end
        if isfield(spec,'an'),  lg{end+1} = '\epsilon'' Analytic'; end
    case 'imag'
        if isfield(spec,'fft'), lg{end+1} = '\epsilon'''' FFT'; end
        if isfield(spec,'an'),  lg{end+1} = '\epsilon'''' Analytic'; end
    case 'cole'
        if isfield(spec,'fft'), lg{end+1} = 'FFT'; end
        if isfield(spec,'an'),  lg{end+1} = 'Analytic'; end
end
end