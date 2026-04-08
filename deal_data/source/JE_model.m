clear
close all
clc
%%
%unit vol and phi is V, dz is nm, charge is e/nm3, dphi is V/nm,
%time and dt is ps
% J unit is e/(nm2 ps)
%ie and q unit is ev/(nm3 ps)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
struct_input = struct();
struct_input.mol = 'PC';
struct_input.T = 298; %K
struct_input.V = 2; %V
struct_input.scanrate = 0; %ps
struct_input.lx = 4.176; %nm
struct_input.ly = 4.254; %nm
struct_input.lz = 10; %nm

analyze_input = struct();
analyze_input.N2mean = 21; %average time
analyze_input.nbin = 1000;
analyze_input.Nmol = 3;
analyze_input.begin_case = 1;
analyze_input.end_case = 50;
analyze_input.dz = struct_input.lz / analyze_input.nbin;

%%%%%%%%%%%%%%%%%%%%%%%%SWITCH%%%%%%%%%%%%%%%%%%%%%%%
global issmooth isfit isforce isAve_onfly isRepeat
issmooth=0;isfit=0;isforce = 0; isAve_onfly=1; isRepeat=1;%repeat control JE
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calc dphi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
analy_mol = struct_input.mol;
a_begin_case = analyze_input.begin_case;
a_end_case = analyze_input.end_case;
onfly_file = sprintf('./onfly_dat/%sonfly%g-%g_dens.dat',analy_mol,a_begin_case,a_end_case);
[charge, charge_d, density, density_d]=remake_onfly_den(onfly_file,struct_input);
disp('success for get onfly data');
den_dt = 0.002*200;...onfly_dt; %ps
skip_t = 500; %ps
analyze_input.dt = den_dt;
onfly_den=struct();
onfly_den.charge=charge;
onfly_den.charge_d=charge_d;    
onfly_den.density=density;
onfly_den.density_d=density_d;
%remake elecharge

Vsdir = sprintf('./CVdat/%s_eleCharge%gV%gps_%g-%g.dat',analy_mol,...
        struct_input.V,struct_input.scanrate,a_begin_case,a_end_case);
Vs=importdata(Vsdir);
Vs_dt = sscanf(Vs.colheaders{1}, '# qout = %f');
Vs_data = Vs.data;
Vs_data = Vs_data/struct_input.lx/struct_input.ly;

elecharge = remake_eleQ(Vs_data,Vs_dt,den_dt,skip_t);
ele_save_dir = sprintf('./CVdat/%seleCharge.mat',struct_input.mol);
save(ele_save_dir,'elecharge');
% t=linspace(0,Vs_dt*length(Vs_data)-skip_t,length(elecharge));
% plot(t,elecharge);hold on
% xlabel('{\itt} (ps)')
% ylabel('{\it-Q_{elec}} (e/nm^2)')
% set(gca,'Units','centimeters','Position',[4 2.5 7 6])
disp('success for get CV data');
%calc dphi
E_dir = sprintf('./Electric_field/%s_elecField%gV%gps.mat',analy_mol,...
        struct_input.V,struct_input.scanrate);
if isforce ~= 1 && exist(E_dir,'file')
    E_data = load(E_dir);
    elec_filed = E_data.elec_filed;
else
    elec_filed = getE(elecharge,onfly_den.charge,struct_input.lz,analyze_input.nbin);
    t=linspace(den_dt,den_dt*size(elec_filed,2),size(elec_filed,2));
    save(E_dir,'elec_filed','t');
end
%calc and save JE
JE_model=getQ(onfly_den, elec_filed , struct_input, analyze_input);

%plot
begin_t = 200;
mdheat_str = sprintf('./mdHeat/%smdHeat%gV%gps_%g-%g.dat',...
    struct_input.mol,struct_input.V,struct_input.scanrate,analyze_input.begin_case,analyze_input.end_case);
mdHeat = xvgread(mdheat_str);
JE_heat = compare_MD(JE_model,mdHeat,analyze_input.dz,struct_input.lx,struct_input.ly);
Idx = (mdHeat(:,1) >= begin_t) & (mdHeat(:,1) <= skip_t);
figure
plot(mdHeat(:,1), mdHeat(:,2), 'Color', [0.8 0.8 0.8]); % 画出整条灰色原线
hold on;
plot(mdHeat(Idx,1), mdHeat(Idx,2), 'r', 'LineWidth', 2); % 用红色高亮你选中的那段

flat_data = mean(mdHeat(Idx,2));
calc_idx = (mdHeat(:,1) > skip_t);
figure
plot(mdHeat(calc_idx,1)-skip_t, -smooth(mdHeat(calc_idx,2),40)+flat_data, 'Color', [0.7 0.7 0.7],'LineWidth',3);
hold on
heat_ion = (intergral(intergral(JEQ_div{2}(1:1000,:),0.01,1,'l'),0.002*200,2)+...
    intergral(intergral(JEQ_div{3}(1:1000,:),0.01,1,'l'),0.002*200,2))*4.176*4.254;
heat_sol = intergral(intergral(JEQ(80:960,:),0.01,1,'l'),0.002*200,2)*4.176*4.254...
    *1.5;
plotIdx = floor(linspace(1,length(JE_heat.time),length(JE_heat.time)/10));
plot(JE_heat.time(plotIdx), heat_ion(plotIdx),'Color',[0 0 0.8],'LineStyle','--')
plot(JE_heat.time(plotIdx), heat_sol(plotIdx),'Color',[0.8 0 0],'LineStyle','--')
plot([0 JE_heat.time(end)],[0 0],'k:','LineWidth',1)
leg=legend('MD','-Q_{all}','-Q_{ion}','Location','north','Orientation','horizontal');
leg.FontSize=14;
     
xlabel('{\itt} (ps)')
ylabel('-{\itQ} (kJ mol^{-1})')
export_fig 'C:\Users\Administrator\Desktop\Ultimate_Figure.png' -png -r600 -m10 -a4 -opengl;

%%
function elecharge = remake_eleQ(Vs_data,Vs_dt,den_dt,skip_t)
global isAve_onfly
    if Vs_dt > den_dt
        error('The time step of electrostatic potential data should be less than that of density data!');
    end
    if mod(den_dt ,Vs_dt) ~= 0
        error('The time step of density data should be integer multiples of that of electrostatic potential data!');
    end
    if(mod(length(Vs_data),10) ~= 0)
        if mod(length(Vs_data)-1,10) ~=0
            error(['The length of the electrode charge is not an integer multiple of 10. Check your CPM_elecharge file ', ...
                 'and ensure it is correct! If there are continuation runs, manual processing will likely be required!',...
                 'If you are confident that you are correct after careful checking, comment out these three lines,',...
                 'and manually modify the code below.']);
        else
            Vs_data = Vs_data(2:end);%去除初始时刻的电荷量
        end
    end
    
    %The ratio of the on-the-fly output frequency to the CPM electrode charge output frequency.
    step_ratio = den_dt / Vs_dt;
    %Skip the non-electrified time period.
    skip_num = skip_t/Vs_dt;
    %Remove step 0 and skip to the step where formal electrification begins.
    Vs_data_skip = Vs_data(skip_num+1 : end);
    %Find the electrode charge at a specific instantaneous step or the average over an interval.
    if isAve_onfly
        num_intervals = floor(length(Vs_data_skip) / step_ratio);
        elecharge = zeros(1, num_intervals);
        for i = 1:num_intervals
            start_idx = (i-1)*step_ratio + 1;
            end_idx = i*step_ratio;
            elecharge(i) = mean(Vs_data_skip(start_idx:end_idx));
        end
    else
        elecharge = Vs_data_skip(step_ratio:step_ratio:length(Vs_data_skip));
    end
end

%
function JE_model=getQ(onfly_den, elec_filed,struct_input, analyze_input)
global issmooth isfit isRepeat
if issmooth==1
    if isfit==1
        save_fold='./Qresult/issmooth/isfit';
    else
        save_fold='./Qresult/issmooth/notfit';
    end
elseif isfit==0
    save_fold='./Qresult/notsmooth/isfit';
else
    save_fold='./Qresult/notsmooth/notfit';
end

savefile=sprintf('%s/%s%gk%gv%gpsIE.mat', save_fold, struct_input.mol,struct_input.T,struct_input.V,struct_input.scanrate);
if isRepeat~=1 && exist(savefile, 'file')
    disp('existing file!');
    JE_model=load(savefile);
    return;
else
    l_E=size(elec_filed,2);
    tmp_elec_filed = ( 0.5*( elec_filed(:,1:l_E-1)+ elec_filed(:,2:l_E) ));
    J=getNrate(onfly_den.charge, analyze_input.N2mean, analyze_input.nbin, analyze_input.dt, analyze_input.dz);
    J_div=getNrate(onfly_den.charge_d, analyze_input.N2mean, analyze_input.nbin, analyze_input.dt, analyze_input.dz);
    
    for i=1:analyze_input.Nmol
        JEQ_div{i}=(-1*J_div{i}.*tmp_elec_filed)*96.488;
    end
    JEQ=(-1*J.*tmp_elec_filed)*96.488;     
end
% if istry==1
%     	tmp=linspace(0.5,1,(1+size(Q,2)-floor(size(Q,2)/2.4)));
%     	tmp=ones(nbin,1)*tmp;
%     	Q(:,floor(size(Q,2)/2.4):size(Q,2))=Q(:,floor(size(Q,2)/2.4):size(Q,2)).*(1-tmp);
% end         %% 2019/3/6 somooth data
disp('complete calculate.save data...');
onfly_dt = analyze_input.dt;
JE_model.JEQ=JEQ;
JE_model.JEQ_div = JEQ_div;
JE_model.J=J;
JE_model.J_div=J_div;
JE_model.onfly_dt=onfly_dt;
save(savefile,'JEQ','JEQ_div','J','J_div','onfly_dt');
end

function Nrate=getNrate(density,ntime,nbin,dt,dz)
    global issmooth isfit
    if ~iscell(density)
         Idensity{1}=density;
    else
        Idensity=density;
    end
    for in=1:length(Idensity)
        if issmooth==1
            Ndensity=meanNsym(Idensity{in},ntime,0);    
        else
            Ndensity=Idensity{in};
        end
        sumdensity2=zeros(nbin,size(Ndensity,2));
        sumdensity=intergral(Ndensity,dz,1);  % for sumdensity%e/nm^2 
        fitnumber=22;
        if isfit==1
            for i=1:nbin                              %fit
                [p,~,mu]=polyfit(1:size(sumdensity,2),sumdensity(i,:),fitnumber); %????????  
                sumdensity2(i,:)=polyval(p,1:size(sumdensity,2),[],mu);
                % p=polyfit(1:size(sumdensity,2),sumdensity(i,:),fitnumber); %????????  
                % sumdensity2(i,:)=polyval(p,1:size(sumdensity,2));
                % sumdensity2(i,:)=smooth(sumdensity(i,:),200,'lowess');
                %   if mod(i,100)==0
                %       i
                %   end
            end
        else
            sumdensity2=sumdensity;
        end
        % plot(sumdensity(40,:)),hold on
        % plot(sumdensity2(40,:));hold off
        %time=0;      
       %%
        diffsumdensity=diff(sumdensity2,1,2);
        Nrate{in}=diffsumdensity/dt*(-1);%%unit is e/nm2/ps 
    end
    if ~iscell(density)
        Nrate=Nrate{1};
    end
end


function [charge, charge_d, density, density_d] = remake_onfly_den(inputFile,struct_input)
global isforce
mol=struct_input.mol;
tem=struct_input.T;
V=struct_input.V;
scanrate=struct_input.scanrate;
savedir=sprintf('./onfly_dat/meandata/%scharge%gk%gV%gps.mat',mol,tem,V,scanrate);

if isforce ~= 1 && exist(savedir,'file')
    den = load(savedir);
    charge = den.charge;
    charge_d = den.charge_d;
    density = den.density;
    density_d = den.density_d;
    return;
end

[ret, ~, ~, ~, ~, ~] = loadOnflyData3D(inputFile);

%Here, the number 4 indicates that onfly will output mass density, number density, charge density and centroid number density.
%The key to successful execution here is that on-the-fly (onfly) calculations are performed in only one region.
for n=1:size(ret{1},2)/4
    %3 and 2 represent charge density and centroid number density, respectively.
    for i=1:length(ret)
        charge_d{n}(:,i) = ret{i}(:,(n-1)*4+3);
        density_d{n}(:,i) = ret{i}(:,(n-1)*4+2);
    end
    if n == 1
        charge = charge_d{1};
        density = density_d{1};
    else
        charge=charge+charge_d{n};
        density=density+density_d{n};
    end
    
end
save(savedir,'charge','charge_d','density','density_d')
end
    
function [ret, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(fnm)
allData = importdata(fnm);
data = allData.data;
textData = allData.textdata;
Lbox = double(regexp(string(textData{2}), "[\d.+-]+", "match"));
lowPos = double(regexp(string(textData{3}), "[\d.+-]+", "match"));
upPos = double(regexp(string(textData{4}), "[\d.+-]+", "match"));
nbin = double(regexp(string(textData{6}), "[\d.+-]+", "match"));
totNbin = double(regexp(string(textData{8}), "[\d.+-]+", "match"));
totFrames = double(regexp(string(textData{9}), "[\d.+-]+", "match"));
xyzRange = cell(1);
for ii=1:length(nbin)
    xyzRange{ii} = linspace(lowPos(ii), upPos(ii), nbin(ii)); 
end
[~, col] = size(data);
nRegion = length(nbin) / 3;
ret = cell(1);
for ii=1:totFrames
    tmpData = data((ii-1)*totNbin+1 : ii*totNbin,:);
    ret{ii} = getRet(tmpData,nRegion,nbin,col);
end
end
    
function ret=getRet(tmpData,nRegion,nbin,col)
ret = cell(1);
tmpIndex = 0;
for ii = 1:nRegion
    eachNbin = nbin(3*(ii-1)+1) * nbin(3*(ii-1)+2) * nbin(3*ii);
    ret{ii} = squeeze(reshape(tmpData((tmpIndex+1):(tmpIndex+eachNbin), :), [nbin((3*(ii-1)+1):3*ii), col]));
    tmpIndex = tmpIndex + eachNbin;
end
if nRegion == 1
    ret = ret{1};
end
end

function dphi=getE(elecharge,charge,lz,nbin)
episi=8.854187817e-21;%F/nm%8.85×10^-12 F/m
dz=lz/nbin;%nm
inter=intergral(charge,dz,1);%e/nm^2
elecharge=elecharge*1.6022e-19;%C/nm^2
inter=inter*1.6022e-19;%C/nm^2
for i=1:length(elecharge)
    dphi(:,i)=-(1/episi)*((inter(:,i)-elecharge(i)));
end
end
% 

%dimension 1 is row, 2 is col
function inter=intergral(fx,dx,dimension,control)
    if dimension==2
        fx=fx';%set the dimension of intergral,x is inter ,y is other 
    end
    [ldx,ldy]=size(fx);
    inter=zeros(ldx,ldy);
    for j=1:ldy
        for i=2:ldx
            inter(i,j)=inter(i-1,j)+(fx(i,j)+fx(i-1,j))*dx/2;
        end
    end
    if exist('control','var')
        inter=inter(size(inter,1),:);
    end
    if dimension==2
        inter=inter';
    end
end
