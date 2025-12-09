clear
close all
clc


dt.water=0.04;
dt.ACN=0.002*50;
dt.PC=0.002*200;
dt.DME=0.002*200;
dt.MT=0.002*50;
dt.AC=0.002*100;
dt.DMSO=0.002*200;

% water=xvgread('water_ave_angleOF_dipole_Z_1-100.dat');
PC=xvgread('PC_ave_angleOF_dipole_Z_1-50.dat');
DMSO=xvgread('DMSO_ave_angleOF_dipole_Z_1-50.dat');
ACN=xvgread('ACN_ave_angleOF_dipole_Z_1-100.dat');
MT=xvgread('MT_ave_angleOF_dipole_Z_1-100.dat');
AC=xvgread('AC_ave_angleOF_dipole_Z_1-100.dat');
DME=xvgread('DME_ave_angleOF_dipole_Z_1-50.dat');
E_water=load('C:\Users\Administrator\Desktop\mycode\JEModel/Electric_field/SOL_elecField2V0ps.mat');
E_PC=load('C:\Users\Administrator\Desktop\mycode\JEModel/Electric_field/PC_elecField2V0ps.mat');
E_DMSO=load('C:\Users\Administrator\Desktop\mycode\JEModel/Electric_field/DMSO_elecField2V0ps.mat');
E_ACN=load('C:\Users\Administrator\Desktop\mycode\JEModel/Electric_field/ACN_elecField2V0ps.mat');
E_MT=load('C:\Users\Administrator\Desktop\mycode\JEModel/Electric_field/MT_elecField2V0ps.mat');
E_AC=load('C:\Users\Administrator\Desktop\mycode\JEModel/Electric_field/AC_elecField2V0ps.mat');
E_DME=load('C:\Users\Administrator\Desktop\mycode\JEModel/Electric_field/DME_elecField2V0ps.mat');
% 创建时间向量
% t_water=linspace(dt.water,dt.water*length(water),length(water));
t_ACN=linspace(dt.ACN,dt.ACN*length(ACN),length(ACN));
t_PC=linspace(dt.PC,dt.PC*length(PC),length(PC));
t_DME=linspace(dt.DME,dt.DME*length(DME),length(DME));
t_MT=linspace(dt.MT,dt.MT*length(MT),length(MT));
t_AC=linspace(dt.AC,dt.AC*length(AC),length(AC));
t_DMSO=linspace(dt.DMSO,dt.DMSO*length(DMSO),length(DMSO));

% 数据可视化 - 绘制所有溶剂的偶极矩角度随时间变化
figure


% plot(t_water, water, 'b-', 'LineWidth', 2);hold on;

plot(t_PC, smooth(PC,20), '-', 'LineWidth', 1.5); hold on;
plot(t_DMSO, smooth(DMSO,20), '-', 'LineWidth', 1.5);
plot(t_ACN, smooth(ACN,20), '-', 'LineWidth', 1.5);
plot(t_MT, smooth(MT,20), '-', 'LineWidth', 1.5);
plot(t_AC, smooth(AC,20), '-', 'LineWidth', 1.5);
plot(t_DME, smooth(DME,20), '-', 'LineWidth', 1.5);
legend('PC','DMSO','ACN','MT','AC','DME')
xlabel('{\itt} (ps)')
ylabel('<\theta_{dipole}> (°)')

figure
plot(E_PC.t,smooth(mean(E_PC.elec_filed(300:700,:)),20))
hold on;
yyaxis right
plot(t_PC, smooth(PC,20), '-', 'LineWidth', 1.5); 
yyaxis left
ylim([0 0.1])
figure
k_pc=smooth(mean(E_PC.elec_filed(300:700,:)),200)./(smooth(PC(2:end),20)-90);
k_water=smooth(mean(E_water.elec_filed(300:700,:)),200)./(smooth(water(2:end),20)-90);
k_dmso=smooth(mean(E_DMSO.elec_filed(300:700,:)),200)./(smooth(DMSO(2:end),20)-90);
k_acn=smooth(mean(E_ACN.elec_filed(300:700,:)),200)./(smooth(ACN(2:end),20)-90);
k_mt=smooth(mean(E_MT.elec_filed(300:700,:)),200)./(smooth(MT(2:end),20)-90);
k_ac=smooth(mean(E_AC.elec_filed(300:700,:)),200)./(smooth(AC(2:end),20)-90);
k_dme=smooth(mean(E_DME.elec_filed(300:700,:)),200)./(smooth(DME(2:end),20)-90);
plot(E_water.t,smooth(k_water,100,'rlowess'));hold on
plot(E_PC.t,smooth(k,100,'rlowess'));
plot(E_DMSO.t,smooth(k_dmso,100,'rlowess'));
plot(E_ACN.t,smooth(k_acn,100,'rlowess'));
plot(E_MT.t,smooth(k_mt,100,'rlowess'));
plot(E_AC.t,smooth(k_ac,100,'rlowess'));
plot(E_DME.t,smooth(k_dme,100,'rlowess'));
legend('water','PC','DMSO','ACN','MT','AC','DME')



