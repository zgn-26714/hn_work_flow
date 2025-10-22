clear
V=2;
scanrate=0;
[ret, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D('/data1/huangnan/ACN/charging/298k/2V/0ps/xvgdata/onfly/onfly1_54.dat_dens.dat');


for n=1:size(ret{1},2)/4
    for i=1:length(ret)
        charge_d{n}(:,i) = ret{i}(:,(n-1)*4+3);
    end
    for i=1:length(ret)
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
plot(charge(:,1))
savedir=sprintf('./meandata/ACNcharge298k%gV%gps_1_54.mat',V,scanrate);
save(savedir,'charge','charge_d','density','density_d')
exit

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
