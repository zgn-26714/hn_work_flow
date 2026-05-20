/* yeting
 * 2020.12.14
 */

/* huangnan
 * modify
 * 2024.11.02
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "include/core_getopt.h"

#ifndef ANALYZE_PRECISION
#define ANALYZE_PRECISION double
#endif
using real = ANALYZE_PRECISION;
using std::vector;
using std::string;

//modify by huangnan,read 3D velocity for every atoms 20240814
struct OnflyPostHandle
{
public:
	OnflyPostHandle( string fnm);
	void readHead();
	bool readNextFrame();

public:
	std::ifstream fid;
	int headSize;
	string commandLine;
	int ngrps;
	int nRegion;
	int userint1;
	int userint2;
	int onfly_flag;
	double dt;
	vector<real> lowPos;
	vector<real> upPos;
	vector<real> dbin;
	vector<int> nbin;
	vector<int> napm;  //add data
	int totNbin;
	real Lbox[3];
	vector<vector<vector<real>>> ivelA, ivelMM, ivelMG, iforceA, iforceMM, iforceMG;
	vector< vector<real>> imassdens, inumdensM, ichargedens, inumdensG;
};

std::string get_inputFile(std::string input) {
	int num, index = -1, len = 0;
	std::string tmp;
	for (int i = 0; i < input.size(); i++) {
		if ('0' <= input[i] && input[i] <= '9') {
			tmp += input[i];
			if (len == 0) index = i;
			len++;
		}
		else {
			if (len != 0) break;
		}
	}
	num = std::stoi(tmp);
	num++;
	if (-1 == index) {
		std::cout << "Error path:\n\t" << input << std::endl;
		exit(-1);
	}
	input.replace(index, len, std::to_string(num));
	return input;
}

vector<vector<vector<vector<real>>>> velA, velMM, velMG, forceA, forceMM, forceMG;
vector<vector<vector<real>>> massdens, numdensM, chargedens, numdensG;
void print_Data(int totNbin, int ngrps, int onfly_flag, int t, std::ofstream& ofile) {
	for (int j = 0; j < totNbin; j++) {
		for (int i = 0; i < ngrps; i++) {
			if ( onfly_flag == 0 ||  onfly_flag == 3 ||  onfly_flag == 4 ||  onfly_flag == 6) {
				ofile << std::setw(15) << massdens[t][i][j] << "\t"
					<< std::setw(15) << numdensM[t][i][j] << "\t"
					<< std::setw(15) << chargedens[t][i][j] << "\t"
					<< std::setw(15) << numdensG[t][i][j] << "\t";
			}
			if ( onfly_flag == 1 ||  onfly_flag == 3 ||  onfly_flag == 5 ||  onfly_flag == 6) {
				ofile << std::setw(15) << velA[t][i][j][0] << "\t" << velA[t][i][j][1] << "\t" << velA[t][i][j][2] << "\t"
					<< std::setw(15) << velMM[t][i][j][0] << "\t" << velMM[t][i][j][1] << "\t" << velMM[t][i][j][2] << "\t"
					<< std::setw(15) << velMG[t][i][j][0] << "\t" << velMG[t][i][j][1] << "\t" << velMG[t][i][j][2] << "\t";
			}
			if ( onfly_flag == 2 ||  onfly_flag == 4 ||  onfly_flag == 5 ||  onfly_flag == 6) {
				ofile << std::setw(15) << forceA[t][i][j][0] << "\t" << forceA[t][i][j][1] << "\t" << forceA[t][i][j][2] << "\t"
					<< std::setw(15) << forceMM[t][i][j][0] << "\t" << forceMM[t][i][j][1] << "\t" << forceMM[t][i][j][2] << "\t"
					<< std::setw(15) << forceMG[t][i][j][0] << "\t" << forceMG[t][i][j][1] << "\t" << forceMG[t][i][j][2] << "\t";
			}
		}
		ofile << '\n';// every bin print '\n'
	}
}




int main(int argc, char** argv)
{
	string inputFnm  = "densMNC3D.onfly";
	string outputFnm = "onfly3D.dat";
	double beginTime = 0;
	double endTime = 0;
	int numFile = 1;
	int dynamic = 1;
    
/*-------------------------get commandline-----------------------------*/
	itp::Getopt getopt(argc, argv, "Onfly PostAnalysis Code");
	getopt.getFixPos(inputFnm, 1, true, "input file name");
	getopt(beginTime, "-b", false, "begin time (ps)");
	getopt(endTime, "-e", false, "end time (ps)");
	getopt(outputFnm, "-o", false, "output file name");
	getopt(numFile, "-n", false, "input file number");
	getopt(dynamic, "-d", false, "isDynamic?");
	getopt.finish();

	std::cout << "\nMake sure you have this stucture:"
		<< "\n\t./case1/file"
		<< "\n\t./case2/file"
		<< "\n\t./case3/file"
		<< "\n\t\t...\n";
	printf("beginTime: %f\n", beginTime);
	printf("endTine: %f\n", endTime);
	printf("input file name: %s\n", inputFnm.c_str());
	printf("output file name: %s\n", outputFnm.c_str());
	printf("input file number: %d\n", numFile);
	printf("isdynamic: %d\n", dynamic);

/*------------------------------initial----------------------------------*/
	OnflyPostHandle ini(inputFnm);
	ini.readHead();

	printf("freq: %d\n", ini.userint2);
	printf("dt: %f\n", ini.dt);

	double dt = ini.dt * ini.userint2;
	size_t beginStep = beginTime / dt;
	size_t endStep = endTime / dt;
	size_t totFrames = 0;

	if (endStep == 0)
	{
		endStep = -1;
	}
	for (size_t i = 0; i < beginStep; i++)
	{
		if (!ini.readNextFrame())
		{
			printf("Error, no more frame!\n");
			std::exit(-1);
		}
	}
	for (size_t i = beginStep; i <= endStep; i++)
	{
		if (!ini.readNextFrame())
			break;
		else
			totFrames++;
	}
	size_t timeindex = totFrames;
	if (0 == dynamic) totFrames = 1;
	std::cout << "totFrames: " <<totFrames << "\n";
	if (endStep == 0)
	{
		endStep = -1;
	}
	
	if (ini.onfly_flag == 0 || ini.onfly_flag == 3 || ini.onfly_flag == 4 || ini.onfly_flag == 6) {
		massdens.resize(totFrames);
		chargedens.resize(totFrames);
		numdensM.resize(totFrames);
		numdensG.resize(totFrames);
	}
	if (ini.onfly_flag == 1 || ini.onfly_flag == 3 || ini.onfly_flag == 5 || ini.onfly_flag == 6) {
		velA.resize(totFrames);
		velMM.resize(totFrames);
		velMG.resize(totFrames);
	}
	if (ini.onfly_flag == 2 || ini.onfly_flag == 4 || ini.onfly_flag == 5 || ini.onfly_flag == 6) {
		forceA.resize(totFrames);
		forceMM.resize(totFrames);
		forceMG.resize(totFrames);
	}

	for (int i = 0; i < totFrames; i++) {
		if (ini.onfly_flag == 0 || ini.onfly_flag == 3 || ini.onfly_flag == 4 || ini.onfly_flag == 6) {
			massdens[i].resize(ini.ngrps);
			chargedens[i].resize(ini.ngrps);
			numdensM[i].resize(ini.ngrps);
			numdensG[i].resize(ini.ngrps);
			for (int j = 0; j < ini.ngrps; j++) {
				massdens[i][j].resize(ini.totNbin, 0.0);
				chargedens[i][j].resize(ini.totNbin, 0.0);
				numdensM[i][j].resize(ini.totNbin, 0.0);
				numdensG[i][j].resize(ini.totNbin, 0.0);
			}
		}
		if (ini.onfly_flag == 1 || ini.onfly_flag == 3 || ini.onfly_flag == 5 || ini.onfly_flag == 6) {
			velA[i].resize(ini.ngrps);
			velMM[i].resize(ini.ngrps);
			velMG[i].resize(ini.ngrps);
			for (int j = 0; j < ini.ngrps; j++)
			{
				velA[i][j].resize(ini.totNbin);
				velMM[i][j].resize(ini.totNbin);
				velMG[i][j].resize(ini.totNbin);
				for (int k = 0; k < ini.totNbin; k++) {
					velA[i][j][k].resize(3, 0.0);
					velMM[i][j][k].resize(3, 0.0);
					velMG[i][j][k].resize(3, 0.0);
				}
			}
		}
		if (ini.onfly_flag == 2 || ini.onfly_flag == 4 || ini.onfly_flag == 5 || ini.onfly_flag == 6) {
			forceA[i].resize(ini.ngrps);
			forceMM[i].resize(ini.ngrps);
			forceMG[i].resize(ini.ngrps);
			for (int j = 0; j < ini.ngrps; j++)
			{
				forceA[i][j].resize(ini.totNbin);
				forceMM[i][j].resize(ini.totNbin);
				forceMG[i][j].resize(ini.totNbin);
				for (int k = 0; k < ini.totNbin; k++) {
					forceA[i][j][k].resize(3, 0.0);
					forceMM[i][j][k].resize(3, 0.0);
					forceMG[i][j][k].resize(3, 0.0);
				}
			}
		}
	}


std::cout<<"initial success!"<<"\n----------\nonfly_flag: "<<ini.onfly_flag<<"\n";
fflush(stdout);

/*----------------calculate from first file-------------------*/
	for (int count = 0; count < numFile; count++) {
		std::cout << inputFnm << std::endl;
		OnflyPostHandle hd(inputFnm);
		hd.readHead();

		for (size_t i = 0; i < beginStep; i++)
		{
			if (!hd.readNextFrame())
			{
				printf("Error, no more frame!\n");
				exit(-1);
			}
		}
		printf("Analysis ...\n");

		//获取帧数据
		for (size_t t = 0; t < timeindex; t++)
		{
			if (!hd.readNextFrame())
				break;
			for (int i = 0; i < hd.ngrps; i++)
			{
				for (int j = 0; j < hd.totNbin; j++)
				{
					//				std::cout<<"where question?\n";
					if (dynamic){
						if (hd.onfly_flag == 0 || hd.onfly_flag == 3 || hd.onfly_flag == 4 || hd.onfly_flag == 6) {
							massdens[t][i][j] += hd.imassdens[i][j];
							chargedens[t][i][j] += hd.ichargedens[i][j];
							numdensM[t][i][j] += hd.inumdensM[i][j];
							numdensG[t][i][j] += hd.inumdensG[i][j];
							//		std::cout<<"no here1\n";
						}
						if (hd.onfly_flag == 1 || hd.onfly_flag == 3 || hd.onfly_flag == 5 || hd.onfly_flag == 6) {
							for (int k = 0; k < 3; k++) {
								velA[t][i][j][k] += hd.ivelA[i][j][k];
								velMM[t][i][j][k] += hd.ivelMM[i][j][k];
								velMG[t][i][j][k] += hd.ivelMG[i][j][k];
							}
							//		std::cout<<"no here2\n";
						}
						if (hd.onfly_flag == 2 || hd.onfly_flag == 4 || hd.onfly_flag == 5 || hd.onfly_flag == 6) {
							for (int k = 0; k < 3; k++) {
								forceA[t][i][j][k] += hd.iforceA[i][j][k];
								forceMM[t][i][j][k] += hd.iforceMM[i][j][k];
								forceMG[t][i][j][k] += hd.iforceMG[i][j][k];
							}
							//		std::cout<<"no here3\n";
						}
					}
					else{
						if (hd.onfly_flag == 0 || hd.onfly_flag == 3 || hd.onfly_flag == 4 || hd.onfly_flag == 6) {
							massdens[0][i][j] += hd.imassdens[i][j];
							chargedens[0][i][j] += hd.ichargedens[i][j];
							numdensM[0][i][j] += hd.inumdensM[i][j];
							numdensG[0][i][j] += hd.inumdensG[i][j];
							//		std::cout<<"no here1\n";
						}
						if (hd.onfly_flag == 1 || hd.onfly_flag == 3 || hd.onfly_flag == 5 || hd.onfly_flag == 6) {
							for (int k = 0; k < 3; k++) {
								velA[0][i][j][k] += hd.ivelA[i][j][k];
								velMM[0][i][j][k] += hd.ivelMM[i][j][k];
								velMG[0][i][j][k] += hd.ivelMG[i][j][k];
							}
							//		std::cout<<"no here2\n";
						}
						if (hd.onfly_flag == 2 || hd.onfly_flag == 4 || hd.onfly_flag == 5 || hd.onfly_flag == 6) {
							for (int k = 0; k < 3; k++) {
								forceA[0][i][j][k] += hd.iforceA[i][j][k];
								forceMM[0][i][j][k] += hd.iforceMM[i][j][k];
								forceMG[0][i][j][k] += hd.iforceMG[i][j][k];
							}
							//		std::cout<<"no here3\n";
						}
					}
					
				}
			}
		}
	
		//average
		if (0 == dynamic) {
			std::cout << "TotFrames: " << timeindex << '\n';
			for (int i = 0; i < hd.ngrps; i++)
			{
				for (int j = 0; j < hd.totNbin; j++)
				{
					if (hd.onfly_flag == 0 || hd.onfly_flag == 3 || hd.onfly_flag == 4 || hd.onfly_flag == 6) {
						massdens[0][i][j] /=  timeindex;
						chargedens[0][i][j] /= timeindex;
						numdensM[0][i][j] /=  timeindex;
						numdensG[0][i][j] /=  timeindex;;
					}
					if (hd.onfly_flag == 1 || hd.onfly_flag == 3 || hd.onfly_flag == 5 || hd.onfly_flag == 6) {
						for (int k = 0; k < 3; k++) {
							velA[0][i][j][k] /=  timeindex;
							velMM[0][i][j][k] /=  timeindex;
							velMG[0][i][j][k] /=  timeindex;
						}
						//		std::cout<<"no here2\n";
					}
					if (hd.onfly_flag == 2 || hd.onfly_flag == 4 || hd.onfly_flag == 5 || hd.onfly_flag == 6) {
						for (int k = 0; k < 3; k++) {
							forceA[0][i][j][k] /=  timeindex;
							forceMM[0][i][j][k] /= timeindex;
							forceMG[0][i][j][k] /= timeindex;
						}
					}
				}
			}
		}
		inputFnm = get_inputFile(inputFnm);
	}

	// average data
	for (int t = 0; t < totFrames; t++) {
		for (int i = 0; i < ini.ngrps; i++) {
			for (int j = 0; j < ini.totNbin; j++) {
				if (ini.onfly_flag == 0 || ini.onfly_flag == 3 || ini.onfly_flag == 4 || ini.onfly_flag == 6) {
					massdens[t][i][j] /= numFile;
					chargedens[t][i][j] /= numFile;
					numdensM[t][i][j] /= numFile;
					numdensG[t][i][j] /= numFile;
					//		std::cout<<"no here1\n";
				}
				if (ini.onfly_flag == 1 || ini.onfly_flag == 3 || ini.onfly_flag == 5 || ini.onfly_flag == 6) {
					for (int k = 0; k < 3; k++) {
						velA[t][i][j][k] /= numFile;
						velMM[t][i][j][k] /= numFile;
						velMG[t][i][j][k] /= numFile;
					}
					//		std::cout<<"no here2\n";
				}
				if (ini.onfly_flag == 2 || ini.onfly_flag == 4 || ini.onfly_flag == 5 || ini.onfly_flag == 6) {
					for (int k = 0; k < 3; k++) {
						forceA[t][i][j][k] /= numFile;
						forceMM[t][i][j][k] /= numFile;
						forceMG[t][i][j][k] /= numFile;
					}
				}
			}
		}
		if (0 == dynamic) break;
	}
		/*-----------------------------get message over-------------------------------------------*/

// printf("test is");
	//输出文件头
	std::ofstream ofile(outputFnm);
	ofile << "% onflyCommandline = \"" << ini.commandLine << "\";\n";
	ofile << "% Lbox = [";
	for (auto&& i : ini.Lbox)
		ofile << i << " ";
	ofile << "];\n";
	ofile << "% lowPos = [";
	for (auto&& i : ini.lowPos)
		ofile << i << " ";
	ofile << "];\n";
	ofile << "% upPos = [";
	for (auto&& i : ini.upPos)
		ofile << i << " ";
	ofile << "];\n";
	ofile << "% fnm = \"" << outputFnm << "\";\n";
	ofile << "% nbin = [";
	for (auto&& i : ini.nbin)
	{
		ofile << i << " ";
	}
	ofile << "];\n";
	ofile << "% onfly_flag = " << ini.onfly_flag << ";\n";
	ofile << "% totNbin = " << ini.totNbin << ";\n";
	ofile << "% totframes = " << totFrames << ";\n";
	ofile << "% ngrps = " << ini.ngrps << ";\n";

	ofile <<R"(% clear
			% close all
			% clc
			% [ret, xyzRange, nbin, lowPos, upPos, Lbox, ngrps, flag] = loadOnflyData3D('./otherdat/onfly3D.dat');
			% 
			%   
			% for i=1:length(ret)
			%     if(mod(size(ret{1}, 2), 3) ~= 0)
			%         for j = 1:ngrps
			%             density_d{j}(:,i) = ret{i}(:,2+(j-1)*4);%2mean numM,1 mass,4 numG
			%             charge_d{j}(:,i) = ret{i}(:,3+(j-1)*4);
			%         end
			%         second_block = ngrps * 4 +1;
			%     else
			%         second_block = 1;
			%     end
			%     
			%     for j = 1:ngrps
			%         for k = 0:2
			%             if (9 * ngrps == size(ret{1}, 2) - second_block + 1)
			%                 if (flag == 1 ||flag == 3)
			%                     vel_d{j, k}(:,i) = ret{i}(:, second_block+k+(j-1)*9);%+3mean MM,+6mean MG
			%                 else
			%                     force_d{j, k}(:,i) = ret{i}(:, second_block+k+(j-1)*9);
			%                 end
			%             else
			%                 vel_d{j, k}(:,i) = ret{i}(:, second_block+k+(j-1)*9);
			%                 force_d{j, k}(:,i) = ret{i}(:, second_block+9*ngrps+k+(j-1)*9);
			%             end
			%         end
			%     end
			%     
			% end
			% 
			% density=zeros(size(density_d{1}));
			% charge=zeros(size(density_d{1}));
			% for k=1:3
			%     vel{k} = zeros(size(vel_d{1,1}));
			%     force{k} = zeros(size(vel_d{1,1}));
			% end
			% 
			% for i = 1:ngrps
			%     if(exist(density_d,'var'))
			%         density = density+density_d{i};
			%         charge = charge + charge_d{i};
			%     end
			%     for k=1:3
			%         if(exist(vel_d,'var'))
			%             vel{k} = vel{k} + vel_d{i,k};
			%             save
			%         end
			%         if(exist(force_d,'var'))
			%             force{k} = force{k} + force_d{i,k};
			%         end
			%     end
			% end
			%     
			% savestr={};
			% if(exist(density_d,'var'))
			%     savestr{end+1} = 'density'; savestr{end+1} = 'density_d';
			%     savestr{end+1} = 'charge'; savestr{end+1} = 'charge_d';
			% end
			% if(exist(vel_d,'var'))
			%    savestr{end+1}='vel'; savestr{end+1} = 'vel_d';
			% end
			% if(exist(force_d,'var'))
			%    savestr{end+1}='force'; savestr{end+1} = 'force_d';
			% end
			% save('./DenVelFor.mat',savestr{:})
			% 
			% function [ret, xyzRange, nbin, lowPos, upPos, Lbox, ngrps, flag] = loadOnflyData3D(fnm)
			% allData = importdata(fnm);
			% data = allData.data;
			% textData = allData.textdata;
			% Lbox = double(regexp(string(textData{2}), "[\d.+-]+", "match"));
			% lowPos = double(regexp(string(textData{3}), "[\d.+-]+", "match"));
			% upPos = double(regexp(string(textData{4}), "[\d.+-]+", "match"));
			% nbin = double(regexp(string(textData{6}), "[\d.+-]+", "match"));
			% flag = double(regexp(string(textData{7}), "[\d.+-]+", "match"));
			% totNbin = double(regexp(string(textData{8}), "[\d.+-]+", "match"));
			% totFrames = double(regexp(string(textData{9}), "[\d.+-]+", "match"));
			% ngrps = double(regexp(string(textData{10}), "[\d.+-]+", "match"));
			% xyzRange = cell(1);
			% for ii=1:length(nbin)
			%    xyzRange{ii} = linspace(lowPos(ii), upPos(ii), nbin(ii)); 
			% end
			% [~, col] = size(data);
			% nRegion = length(nbin) / 3;
			% ret = cell(1);
			% for ii=1:totFrames
			%     tmpData = data((ii-1)*totNbin+1 : ii*totNbin,:);
			%     ret{ii} = getRet(tmpData,nRegion,nbin,col);
			% end
			% end
			% 
			% function ret=getRet(tmpData,nRegion,nbin,col)
			% ret = cell(1);
			% tmpIndex = 0;
			% for ii = 1:nRegion
			%     eachNbin = nbin(3*(ii-1)+1) * nbin(3*(ii-1)+2) * nbin(3*ii);
			%     ret{ii} = squeeze(reshape(tmpData((tmpIndex+1):(tmpIndex+eachNbin), :), [nbin((3*(ii-1)+1):3*ii), col]));
			%     tmpIndex = tmpIndex + eachNbin;
			% end
			% if nRegion == 1
			%     ret = ret{1};
			% end
			% end)";
	ofile << "\n";
	//ofile.setf(std::ios::scientific);
	//ofile.precision(8);
	
	if (!dynamic) {
		print_Data( ini.totNbin, ini.ngrps, ini.onfly_flag ,0 , ofile);
	}
	else
	{
		for (size_t t = 0; t < totFrames; t++) {
			print_Data(ini.totNbin, ini.ngrps, ini.onfly_flag, t, ofile);
		}
	}
	ofile.close();
	printf("Analysis done.\n");
}

OnflyPostHandle::OnflyPostHandle( string fnm)
{
	fid.open(fnm,  std::ios::binary);
	if (!fid)
	{
		printf("Cannot open file \"%s\", exit!\n", fnm.c_str());
		 exit(-1);
	}
}

void OnflyPostHandle::readHead()
{
	int magicNumber = 0;
	fid.read((char*)&magicNumber, sizeof(int));
	if (magicNumber == 20080513)
	{//please Make sure this corresponds to the other one
        int commandLineSize = 0;
		fid.read((char*)&headSize, sizeof(int));
		fid.read((char*)&commandLineSize, sizeof(int));
		fid.read((char*)&ngrps, sizeof(int));
		fid.read((char*)&nRegion, sizeof(int));
		fid.read((char*)&userint1, sizeof(int));
		fid.read((char*)&userint2, sizeof(int));
		fid.read((char*)&onfly_flag, sizeof(int));
		fid.read((char*)&dt, sizeof(double));
		lowPos.resize(3 * nRegion);
		upPos.resize(3 * nRegion);
		nbin.resize(3 * nRegion);
		dbin.resize(3 * nRegion);
        napm.resize(ngrps);
		fid.read((char*)lowPos.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)upPos.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)dbin.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)nbin.data(), sizeof(int) * 3 * nRegion);
		fid.read((char*)napm.data(), sizeof(int) * ngrps);//add read
		fid.read((char*)Lbox, sizeof(real) * 3);
        commandLine.resize(commandLineSize);
		fid.read((char*)commandLine.data(), commandLine.size());
	}
	else
	{
		 std::cout << "Error format, cannot read head!" <<  std::endl;
		 exit(-1);
	}
	totNbin = 0;
	for (int i = 0; i < nRegion; i++)
		totNbin += nbin[i * 3] * nbin[i * 3 + 1] * nbin[i * 3 + 2];

	
	if (onfly_flag == 0 || onfly_flag == 3 || onfly_flag == 4 || onfly_flag == 6) {
		imassdens.resize(ngrps);
		ichargedens.resize(ngrps);
		inumdensM.resize(ngrps);
		inumdensG.resize(ngrps);
	}
	if (onfly_flag == 1 || onfly_flag == 3 || onfly_flag == 5 || onfly_flag == 6) {
		ivelA.resize(ngrps);
		ivelMM.resize(ngrps);
		ivelMG.resize(ngrps);
	}
	if (onfly_flag == 2 || onfly_flag == 4 || onfly_flag == 5 || onfly_flag == 6) {
		iforceA.resize(ngrps);
		iforceMM.resize(ngrps);
		iforceMG.resize(ngrps);
	}

	for (int i = 0; i < ngrps; i++) {
		if (onfly_flag == 0 || onfly_flag == 3 || onfly_flag == 4 || onfly_flag == 6) {
			imassdens[i].resize(totNbin);
			ichargedens[i].resize(totNbin);
			inumdensM[i].resize(totNbin);
			inumdensG[i].resize(totNbin);
		}
		if (onfly_flag == 1 || onfly_flag == 3 || onfly_flag == 5 || onfly_flag == 6) {
			ivelA[i].resize(totNbin);
			ivelMM[i].resize(totNbin);
			ivelMG[i].resize(totNbin);
			for (int j = 0; j < totNbin; j++)
			{
				ivelA[i][j].resize(3);
				ivelMM[i][j].resize(3);
				ivelMG[i][j].resize(3);
			}
		}
		if (onfly_flag == 2 || onfly_flag == 4 || onfly_flag == 5 || onfly_flag == 6) {
			iforceA[i].resize(totNbin);
			iforceMM[i].resize(totNbin);
			iforceMG[i].resize(totNbin);
			for (int j = 0; j < totNbin; j++)
			{
				iforceA[i][j].resize(3);
				iforceMM[i][j].resize(3);
				iforceMG[i][j].resize(3);
			}
		}
	}
}

bool OnflyPostHandle::readNextFrame()
{
	int magicNumber = 0;
	int succeess = 0;
	while (true)
	{
		fid.read((char*)&magicNumber, sizeof(int));
		if (magicNumber == 20241101)
		{
			for (int i = 0; i < ngrps; i++)
			{
				fid.read((char*)(imassdens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(inumdensM[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(ichargedens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(inumdensG[i].data()), sizeof(real) * totNbin);
			}
			succeess = 1;
		}
		else if (magicNumber == 20080513)
		{
			int ignoreSize = 0;
			fid.read((char*)&ignoreSize, sizeof(int));
			fid.seekg(ignoreSize - 2 * sizeof(int), std::ios::cur);
			magicNumber = 0;
			continue;
		}
		else if (magicNumber == 20241102)
		{
			for (int i = 0; i < ngrps; i++)
			{
				for (int j = 0; j < totNbin; j++) {
					fid.read((char*)(ivelA[i][j].data()), sizeof(real) * 3);
					fid.read((char*)(ivelMM[i][j].data()), sizeof(real) * 3);
					fid.read((char*)(ivelMG[i][j].data()), sizeof(real) * 3);
				}
			}
			succeess = 1;
		}
		else if (magicNumber == 20241103) {
			for (int i = 0; i < ngrps; i++)
			{
				for (int j = 0; j < totNbin; j++) {
					fid.read((char*)(iforceA[i][j].data()), sizeof(real) * 3);
					fid.read((char*)(iforceMM[i][j].data()), sizeof(real) * 3);
					fid.read((char*)(iforceMG[i][j].data()), sizeof(real) * 3);
				}
			}
			succeess = 1;
		}
		else if (magicNumber == 20240000 && succeess == 1) {
			return true;
		}
		else {
			return false;
		}
	}
}