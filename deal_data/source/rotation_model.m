clc
clear
close all
%% parameters
p=0.0489; % unit e nm?diplo of water molecular
p = -0.088593706414768; %e nm,dipole of ACN
KbT=0.02568; % unit is eV?thermal energy
% zeta=0.01;% unit is eV ps 
I = 0.0194 ; % Moment of Inertia in amu nm^2  (0.0194)
I = 0.5674;  %ACN
I=I *1.6606390666e-27*1e-18;% Moment of Inertia in kg m^2
I = I * 1e24 * 6.2415e18; % Convert kg m^2 to eV ps^2
Elistequ=linspace(0,0.2,100);% Electric filed, unit is V/nm
%% Equiliburim E-costheta
% getEqubrm=@getThetaLgV;
getEqubrm=@getEqubrimMD;
ThetaEqu=getEqubrm(Elistequ,p,KbT);
% zetas=[0.01 0.1 0.2 0.4 0.8 1.6];% unit is eV ps 闃诲凹绯绘暟
zetas=[0.23];% unit is eV ps 闃诲凹绯绘暟
%% get E(t)
[Et,tlistEt]=getEt;
%%
% initial condition
theta0=pi/2;
omega0 = 0;
Y0 = [theta0; omega0];
tspan = [0 tlistEt(end)];
figure;
for iz=1:length(zetas)
    zeta=zetas(iz);
    ODEsolve=@(t,Y)myODE(t, Y,I,zeta,p,Et,tlistEt,ThetaEqu,Elistequ);
    % resolve
    [t, Y] = ode45(ODEsolve, tspan, Y0);
    thetasolve=Y(:,1);
    thetasolve=thetasolve/pi*180;% rad-> degree
    omegasolve=Y(:,2);
    omegasolve=omegasolve/pi*180;% rad/2-> degree/2
    plot(t, thetasolve, '-'),hold on
    xlabel('t (ps)');ylabel('\theta (^o)');
    xlim([0 200])
    legendlabels{iz}=['\zeta=',num2str(zeta)];
    tmp=diff(thetasolve/180*pi)./diff(t);
    current=-tmp*0.0489*32.8;
    save('SOLcurrent.mat','current','t')
    save([num2str(zeta),'data_2.mat'],'t','thetasolve')
end
legend(legendlabels,'fontsize',16)


function dYdt = myODE(t, Y,I,zeta,p,Et,tlistEt,ThetaEqu,Eequ)
theta = Y(1);  % 绗竴涓彉閲忔槸胃
omega = Y(2);  % 绗簩涓彉閲忔槸蠅
Enowequ=interp1(ThetaEqu,Eequ,theta,'linear','extrap');
if isnan(Enowequ) || theta<min(ThetaEqu)  ||  theta>max(ThetaEqu)
    warning('interp may wrong')
end
Enowt=interp1(tlistEt,Et,t);
if isnan(Enowt)  || theta<min(ThetaEqu)  ||  theta>max(ThetaEqu)
    warning('interp may wrong','linear','extrap')
end
% dot_omega = (-zeta * omega + (Enowt-Enowequ) * p * sin(theta)) / I;
dot_omega = (-zeta * omega + (Enowt-Enowequ) * p * pi/4) / I;%test
dYdt = [omega; dot_omega];

end




function [Et,tlistEt]=getEt
dt=0.04;
dat=load('./MDdata/gap10dphi298k2V0ps.mat');% E unit is V/nm ; theta unit is dgree
Emesh=dat.dphi;
[Nz,Nt]=size(Emesh);
sz=floor(Nz*0.5/5);
ez=floor(Nz*4.5/5);
Et=mean(Emesh(sz:ez,:));% unit is V/nm
Et(1)=0.195;
Et=[0.2;Et(:)];
tlistEt=linspace(0.04,Nt*dt,Nt); % unit is ps
tlistEt=[0;tlistEt(:)];
end



function theta=getThetaLgV(Elist,p,KbT)
cosTheta=coth(p.*Elist./KbT)-KbT./(p.*Elist);
theta=acos(cosTheta);
end

% function theta=getEqubrimMD(Elist,p,KbT)
% dat=load('./MDdata/Ethetanacl.mat');% E unit is V/nm ; theta unit is dgree
% EMD=[0;dat.EFs(:)];
% MD_cos=[0;dat.MD_cos(:)];
% cosTheta=interp1(EMD,MD_cos,Elist);
% theta=acos(cosTheta);
% theta(1)=pi/2;
% end

function theta=getEqubrimMD(Elist,p,KbT)
dat=load('./MDdata/ACN_Ethetanacl.mat');% E unit is V/nm ; theta unit is dgree
EMD=[dat.EFs(:)];
% MD_cos=[0;dat.MD_cos(:)];
MD_theta=[dat.angle(:);90]/180*pi;
% cosTheta=interp1(EMD,MD_cos,Elist);
theta=interp1(EMD,MD_theta,Elist);
% theta=acos(cosTheta);
% thetaMD=acos(MD_cos);
theta(1)=pi/2;
end%和z轴组成的平面的法线绕着(看我的getI代码)