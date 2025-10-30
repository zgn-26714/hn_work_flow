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
struct_input.mol = getenv('analyze_mol');;
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
analyze_input.dz = struct_input.lz / analyze_input.nbin;

global issmooth isfit %switch
issmooth=1;isfit=1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calc dphi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
onfly_file = getenv('onfly_dat');
[charge, charge_d, density, density_d]=remake_onfly_den(onfly_file);
onfly_den=struct();
onfly_den.charge=charge;
onfly_den.charge_d=charge_d;    
onfly_den.density=density;
onfly_den.density_d=density_d;
%remake elecharge
analy_mol = getenv('analyze_mol');
a_begin_case = getenv('analyze_begin_case');
a_end_case = getenv('analyze_end_case');
Vs=importdata(sprintf('./deal_data/eleQ/%s_eleCharge2V0ps_%s-%s.dat',analy_mol,a_begin_case,a_end_case));
Vs_dt = Vs.data;
Vs_data = str2double(Vs.textdata(2:end, 1));
Vs_data = Vs_data/struct_input.lx/struct_input.ly;
den_dt = str2double(getenv('den_dt')); %ps
skip_t = 400; %ps
elecharge = remake_eleQ(Vs_data,Vs_dt,den_dt,skip_t);

%calc dphi
elec_filed = getE(elecharge,onfly_den.charge,struct_input.lz,analyze_input.nbin);

%calc and save JE
getQ(onfly_den, elec_filed , struct_input, analyze_input);

%%
function elecharge = remake_eleQ(Vs_data,Vs_dt,den_dt,skip_t)
    if Vs_dt > den_dt
        error('The time step of electrostatic potential data should be less than that of density data!');
    end
    if mod(den_dt ,Vs_dt) ~= 0
        error('The time step of density data should be integer multiples of that of electrostatic potential data!');
    end

    step_ratio = den_dt / Vs_dt;
    skip_num = skip_t/Vs_dt;
    Vs_data_skip = Vs_data(skip_num+1 : end);
    elecharge = Vs_data_skip(1:step_ratio:end);
end

%
function getQ(onfly_den, elec_filed,struct_input, analyze_input)
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
if exist(savefile, 'file')
    disp('existing file!');
    return;
else
    l_E=size(elec_filed,2);
    tmp_elec_filed = ( 0.5*( elec_filed(:,1:l_E-1)+ elec_filed(:,2:l_E) ));
    J=getNrate(onfly_den.charge, analyze_input.ntime, analyze_input.nbin, analyze_input.dt, analyze_input.dz);
    J_div=getNrate(onfly_den.charge_d, analyze_input.ntime, analyze_input.nbin, analyze_input.dt, analyze_input.dz);
    for i=1:Nmol
        JEQ_div{i}=(-1*J_div{i}.*tmp_elec_filed)*96.488;
    end
    JEQ=(-1*J.*dphi)*96.488;     
end
% if istry==1
%     	tmp=linspace(0.5,1,(1+size(Q,2)-floor(size(Q,2)/2.4)));
%     	tmp=ones(nbin,1)*tmp;
%     	Q(:,floor(size(Q,2)/2.4):size(Q,2))=Q(:,floor(size(Q,2)/2.4):size(Q,2)).*(1-tmp);
% end         %% 2019/3/6 somooth data
save(savefile,'JEQ','JEQ_div','J','J_div','analyze_input.dt');
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
        fitnumber=20;
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


function [charge, charge_d, density, density_d] = remake_onfly_den(inputFile)
[ret, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(inputFile);

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

mol=getenv('analy_mol');
tem=getenv('analyze_T');
V=str2double(getenv('analyze_V'));
scanrate=str2double(getenv('analyze_tau'));
savedir=sprintf('./deal_data/onfly/meandata/%scharge%sk%gV%gps.mat',mol,tem,V,scanrate);
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
episi=8.85e-21;%F/nm%8.85×10^-12 F/m
dz=lz/nbin;%nm
inter=intergral(charge,dz,1);%e/nm^2
elecharge=elecharge*1.6022e-19;%C/nm^2
inter=inter*1.6022e-19;%C/nm^2
for i=1:size(charge,2)
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
