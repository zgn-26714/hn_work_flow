clc
clear
close all
%% parameters
p=struct();
p.water=0.0489; % unit e nm?diplo of water molecular
p.DMSO=0.103404;
p.ACN=0.082529;
p.PC=0.118873;
p.DME=0.000501;
p.MT=0.052467;
p.AC=0.061040;
% zeta=0.01;% unit is eV ps
I=struct();
% I1 = 0.005881 e amu nm^2
% I2 = 0.013527 e amu nm^2
% I3 = 0.019409 e amu nm^2
I.water = 0.0194 ; % Moment of Inertia in amu nm^2  (0.0194)
% I1 = 0.729713 amu nm^2
% I2 = 0.738537 amu nm^2
% I3 = 1.218290 amu nm^2
I.DMSO=1.218290;
% I1 = 0.045162 amu nm^2
% I2 = 0.527158 amu nm^2
% I3 = 0.536481 amu nm^2
I.ACN=0.536481;
% I1 = 0.983173 amu nm^2
% I2 = 2.048176 amu nm^2
% I3 = 2.500885 amu nm^2
I.PC = 2.500885;
% I1 = 0.269635 amu nm^2
% I2 = 4.004545 amu nm^2
% I3 = 4.145344 amu nm^2
I.DME=4.145344;
% I1 = 0.040079 amu nm^2
% I2 = 0.201462 amu nm^2
% I3 = 0.209210 amu nm^2
I.MT=0.209210;
% I1 = 0.500478 amu nm^2
% I2 = 1.027096 amu nm^2
% I3 = 0.591936 amu nm^2
I.AC=1.027096;

dt = struct();
dt.water=0.04;
dt.ACN=0.002*50;
dt.PC=0.002*200;
dt.DME=0.002*200;
dt.MT=0.002*50;
dt.AC=0.002*100;
dt.DMSO=0.002*200;
%重新写一下相乘（根据自己选的内容：ACN、water？）
name = 'DME';
%AC:0.23/0.25
%MT:0.34
%ACN:0.32 0.33 0.34
%DMSO:~3.75 (3.96 3.56)
%PC:~3.76 (3.59 3.96)
%DME:
p = p.(name);
I = I.(name);
dt=dt.(name);
I = I *1.6606390666e-27*1e-18;% Moment of Inertia in kg m^2
I = I * 1e24 * 6.2415e18; % Convert kg m^2 to eV ps^2
KbT=0.02568; % unit is eV?thermal energy
Elistequ=linspace(0,0.2,100);% Electric filed, unit is V/nm
%% Equiliburim E-costheta
% getEqubrm=@getThetaLgV;
getEqubrm=@getEqubrimMD;
ThetaEqu=getEqubrm(Elistequ,name);
% zetas=[0.01 0.1 0.2 0.4 0.8 1.6];% unit is eV ps 闃诲凹绯绘暟
zetas = [0.003 0.004];% unit is eV ps 闃诲凹绯绘暟
zeta_step=0.0001;
%% get E(t)
[Et,tlistEt]=getEt(name,dt);
%%
% initial condition
theta0=pi/2;
omega0 = 0;
Y0 = theta0;
% Y0 = [theta0; omega0];
tspan = [0 tlistEt(end)];
% figure;
count=1;
MD_theta=xvgread(sprintf('./MDtheta/%s_ave_angleOF_dipole_Z_1-50.dat',name));
MD_t=linspace(0,dt*(length(MD_theta)-1),length(MD_theta));
figure
plot(MD_t,MD_theta,'Color',[0.5 0.5 0.5]);hold on
for zeta= (zetas(1) :zeta_step  : zetas(2))
%     ODEsolve=@(t,Y)myODE(t, Y,I, zeta,p,Et,tlistEt,ThetaEqu,Elistequ);
    ODEsolve=@(t,Y)myODE(t, Y, zeta,p,Et,tlistEt,ThetaEqu,Elistequ);
    % resolve
    options = odeset('MaxStep',0.5) ;
    [t, Y] = ode45(ODEsolve, tspan, Y0, options);
    thetasolve=Y(:,1);
    thetasolve=thetasolve/pi*180;% rad-> degree
%     omegasolve=Y(:,2);
%     omegasolve=omegasolve/pi*180;% rad/2-> degree/2
    all_t(:,count)=t;
    all_thetas(:,count)=thetasolve;
    [MAE(count), MAPE(count), RMSE(count),R2(count)] = calc_error(t,thetasolve, MD_t, MD_theta);
    plot(t, thetasolve,'LineWidth',0.5),hold on
    xlabel('t (ps)');ylabel('\theta (^o)');
    xlim([0 200])
    legendlabels{count}=['\zeta=',num2str(zeta)];
    count=count+1;
%     p=polyfit(1:length(thetasolve),thetasolve',20); %????????  
%     thetasolve=polyval(p,1:length(thetasolve));
%     tmp=diff(thetasolve/180*pi)./diff(t');
%     current=-tmp*0.0489*32.8*8/(3*pi);
    save('0.263data.mat','t','thetasolve')
%     t=t(1:end-1)';
%     save('SOLcurrent.mat','current','t')
end
legendlabels = [{'MD_theta'}, legendlabels];
legend(legendlabels,'fontsize',14)
figure
xindex = (zetas(1) :zeta_step  : zetas(2));
plot(xindex,  MAE/min(MAE));hold on
plot(xindex,  MAPE/min(MAPE));
plot(xindex, RMSE/min(RMSE));
plot(xindex, R2)
[~,b]=max(R2);
plot([xindex(b) xindex(b)],[0.98 1.1],'k--','LineWidth',1);hold on
[~,b]=min(MAE);
plot([xindex(b) xindex(b)],[0.98 1.1],'k--','LineWidth',1);hold on
[~,b]=min(MAPE);
plot([xindex(b) xindex(b)],[0.98 1.1],'k--','LineWidth',1);hold on
[~,b]=min(RMSE);
plot([xindex(b) xindex(b)],[0.98 1.1],'k--','LineWidth',1);hold on
ylim([0.98 1.06])
legend('MAE','MAPE','RMSE','R2')
figure
tmp=abs(xindex-b);
[~,b]=min(tmp);
plot(MD_t,MD_theta,'Color',[0.5 0.5 0.5]);hold on
plot(all_t(:,b),all_thetas(:,b),'k--')
xlabel('t (ps)');ylabel('\theta (^o)');
xlim([0 500])
legend('MD','theory')
set(gca,'Units','centimeters','Position',[4 2.5 7 6])

function dYdt = myODE(t_bar, Y_bar,zeta,p,Et,tlistEt,ThetaEqu,Eequ)
theta = Y_bar;
Enowequ0=interp1(ThetaEqu,Eequ,theta,'linear','extrap');
% if isnan(Enowequ0) || theta<min(ThetaEqu)  ||  theta>max(ThetaEqu)
if isnan(Enowequ0) ||  theta>max(ThetaEqu)
    warning('interp may wrong1')
    theta
end
Enowt0=interp1(tlistEt,Et,t_bar);
if isnan(Enowt0)  ||  theta>max(ThetaEqu)
    warning('interp may wrong2','linear','extrap')
end
omega =  ((Enowt0-Enowequ0) *p* pi/4)/zeta;
% omega =  ((Enowt0-Enowequ0) *p*  sin(theta))/zeta;
dYdt = omega;
if mod(t_bar,1e10)==0
    t_bar
end
end

function [MAE,MAPE,RMSE,R2]= calc_error(t,thetasolve, MD_t, MD_theta)
    % 将模拟角度插值到实验时间点上
    theta_sim_interp = interp1(t, thetasolve, MD_t, 'linear', 'extrap');
    
    % 计算 err
    MAE = sum(abs(theta_sim_interp - MD_theta'));
    MAPE = sum(abs((theta_sim_interp - MD_theta') ./ MD_theta')) * 100;
    RMSE = sqrt(sum((theta_sim_interp - MD_theta').^2));
    
    % R^2
    SS_res = sum((theta_sim_interp - MD_theta').^2);
    SS_tot = sum((theta_sim_interp - (mean(MD_theta))').^2);
    R2 = 1 - SS_res / SS_tot;

end





% function dYdt = myODE(t, Y,I,zeta,p,Et,tlistEt,ThetaEqu,Eequ)
% 
% theta = Y(1);      % 角度 θ
% omega = Y(2);      % 角速度 ω
% 
% %% ---- 计算 E_eq(theta) ----
% % 正确插值：用 θ 查表得到 E_eq
% Enowequ = interp1(ThetaEqu, Eequ, theta, 'linear', 'extrap');
% 
% if theta < min(ThetaEqu) || theta > max(ThetaEqu)
%     warning('theta out of equilibrium table range.');
% end
% 
% %% ---- 计算 E(t) ----
% Enowt = interp1(tlistEt, Et, t, 'linear');
% 
% %% ---- 使用你的公式 ----
% % I * dω/dt = -ζω + (E(t) - Eeq(θ)) * p * π/4
% dot_omega = (-zeta * omega + (Enowt - Enowequ) * p * (pi/4)) / I;
% 
% %% ---- output ----
% dYdt = [omega; dot_omega];
% 
% end



function [Et,tlistEt]=getEt(name,dt)
dat=load(sprintf('./MDdata/%s_elecField2V0ps.mat',name));% E unit is V/nm ; theta unit is dgree
Emesh=dat.elec_filed;
[Nz,Nt]=size(Emesh);
sz=floor(Nz*1/5);
ez=floor(Nz*4/5);
Et=mean(Emesh(sz:ez,:));% unit is V/nm
% Et(1)=0.195;
Et=[0.2;Et(:)];
tlistEt=linspace(dt,Nt*dt,Nt); % unit is ps
tlistEt=[0;tlistEt(:)];
end



function theta=getThetaLgV(Elist,p,KbT)
cosTheta=coth(p.*Elist./KbT)-KbT./(p.*Elist);
theta=acos(cosTheta);
end
% 
% function theta=getEqubrimMD(Elist,p,KbT)
% dat=load('./MDdata/Ethetanacl.mat');% E unit is V/nm ; theta unit is dgree
% EMD=[0;dat.EFs(:)];
% MD_cos=[0;dat.MD_cos(:)];
% cosTheta=interp1(EMD,MD_cos,Elist);
% theta=acos(cosTheta);
% theta(1)=pi/2;
% end

function theta=getEqubrimMD(Elist,name)
dat=load(sprintf('./MDdata/%s_Ethetanacl.mat',name));% E unit is V/nm ; theta unit is dgree
EMD=dat.EFs(:);
% MD_cos=[0;dat.MD_cos(:)];
MD_theta=dat.angle(:)/180*pi;
% cosTheta=interp1(EMD,MD_cos,Elist);
theta=interp1(EMD,MD_theta,Elist);
% theta=acos(cosTheta);
% thetaMD=acos(MD_cos);
theta(1)=pi/2;
end