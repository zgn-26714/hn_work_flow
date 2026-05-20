/* yeting
 * 2020.12.14
 */

//huangnan modified

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "include/core_getopt.h"
#ifndef ANALYZE_PRECISION
#define ANALYZE_PRECISION float
#endif
using real = ANALYZE_PRECISION;



struct OnflyPostHandle
{
public:
	OnflyPostHandle(std::string fnm);
	void readHead();
	bool readNextFrame();

public:
	std::ifstream fid;
	int headSize;
	std::string commandLine;
	int ngrps;
	int nRegion;
	int userint1;
	int userint2;
	double dt;
	std::vector<real> lowPos;
	std::vector<real> upPos;
	std::vector<real> dbin;
	std::vector<int> nbin;
	int totNbin;
	real Lbox[3];
	std::vector<std::vector<real>> massdens;
	std::vector<std::vector<real>> numdensM;
	std::vector<std::vector<real>> chargedens;
	std::vector<std::vector<real>> numdensG;
};

std::string get_inputFile(std::string input) {
	int num, index = -1, len=0;
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


int main(int argc, char** argv)
{
	std::string inputFnm = "densMNC3D.onfly";
	std::string outputFnm = "onfly3D.dat";
	double beginTime = 0;
	double endTime = 0;
	double numFile = 1;
	int dynamic = 1;

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
		<< "\n\t\t...";
	printf("beginTime: %f\n", beginTime);
	printf("endTine: %f\n", endTime);
	printf("input file name: %s\n", inputFnm.c_str());
	printf("output file name: %s\n", outputFnm.c_str());
	printf("input file number: %f\n", numFile);
	printf("isdynamic: %d\n", dynamic);

	OnflyPostHandle ini(inputFnm);
	ini.readHead();

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
	if (0 == dynamic) totFrames = 1;
	std::cout << totFrames << "\n";
	std::vector<std::vector<std::vector<real>>> massdens;
	std::vector<std::vector<std::vector<real>>> numdensM;
	std::vector<std::vector<std::vector<real>>> chargedens;
	std::vector<std::vector<std::vector<real>>> numdensG;
	massdens.resize(totFrames);
	chargedens.resize(totFrames);
	numdensM.resize(totFrames);
	numdensG.resize(totFrames);
	for (int i = 0; i < totFrames; i++){
		massdens[i].resize(ini.ngrps);
		chargedens[i].resize(ini.ngrps);
		numdensM[i].resize(ini.ngrps);
		numdensG[i].resize(ini.ngrps);
		for (int j = 0; j < ini.ngrps; j++){
			massdens[i][j].resize(ini.totNbin, 0);
			numdensM[i][j].resize(ini.totNbin, 0);
			chargedens[i][j].resize(ini.totNbin, 0);
			numdensG[i][j].resize(ini.totNbin, 0);
		}
	}//initial



	for (int count = 0; count < numFile; count++) {
		std::cout << inputFnm << std::endl;
		OnflyPostHandle hd(inputFnm);
		hd.readHead();

		for (size_t i = 0; i < beginStep; i++)
		{
			if (!hd.readNextFrame())
			{
				printf("Error, no more frame!\n");
				std::exit(-1);
			}
		}

		printf("Analysis ...\n");
		size_t nowTime = 0;
		for (size_t t = 0; t < totFrames; t++)
		{
			if (!hd.readNextFrame())
				break;
			#pragma omp parallel for num_threads(8)
			for (int i = 0; i < hd.ngrps; i++)
			{
				for (int j = 0; j < hd.totNbin; j++)
				{
					massdens[t][i][j] += hd.massdens[i][j];
					numdensM[t][i][j] += hd.numdensM[i][j];
					chargedens[t][i][j] += hd.chargedens[i][j];
					numdensG[t][i][j] += hd.numdensG[i][j];
				}
			}
			//std::cout << nowTime << ": \n"; fflush(stdout);
			nowTime++;
			//if( int(nowTime * dt)  % 100 == 0)
			//printf("Read frame from %.3f ps, %zd frames are handled\n",(beginStep + nowTime) * dt, nowTime);
		}
		if (dynamic == 0) {
			std::cout << "TotFrames: " << nowTime << '\n';
			for (int i = 0; i < hd.ngrps; i++)
			{
				for (int j = 0; j < hd.totNbin; j++)
				{
					massdens[0][i][j] /= nowTime;
					numdensM[0][i][j] /= nowTime;
					chargedens[0][i][j] /= nowTime;
					numdensG[0][i][j] /= nowTime;
				}
			}
		}
		inputFnm = get_inputFile(inputFnm);
	}
	
	#pragma omp parallel for num_threads(8)
	for (size_t t = 0; t < totFrames; t++)
	{
		for (int i = 0; i < ini.ngrps; i++)
		{
			for (int j = 0; j < ini.totNbin; j++)
			{
				massdens[t][i][j] /= numFile;
				numdensM[t][i][j] /= numFile;
				chargedens[t][i][j] /= numFile;
				numdensG[t][i][j] /= numFile;
			}
		}
	}

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
	ofile << "% nRegion = " << ini.nRegion << ";\n";
	ofile << "% totNbin = " << ini.totNbin << ";\n";
	ofile << "% totframes = " << totFrames << ";\n";
	ofile << "% [data, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(fnm);\n";
	ofile << "%\n";
	ofile <<
		R"(% function [ret, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(fnm)
% allData = importdata(fnm);
% data = allData.data;
% textData = allData.textdata;
% Lbox = double(regexp(string(textData{2}), "[\d.+-]+", "match"));
% lowPos = double(regexp(string(textData{3}), "[\d.+-]+", "match"));
% upPos = double(regexp(string(textData{4}), "[\d.+-]+", "match"));
% nbin = double(regexp(string(textData{6}), "[\d.+-]+", "match"));
% totNbin = double(regexp(string(textData{8}), "[\d.+-]+", "match"));
% totFrames = double(regexp(string(textData{9}), "[\d.+-]+", "match"));
% xyzRange = cell(1);
% for ii=1:length(nbin)
%    xyzRange{ii} = linspace(lowPos(ii), upPos(ii), nbin(ii)); 
% end
% [~, col] = size(data);
% nRegion = length(nbin) / 3;
% ret = cell(1);
% for ii=1:totFrames
%     tmpData = data((ii-1)*totNbin : ii*totNbin,:);
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
	ofile.setf(std::ios::scientific);
	ofile.precision(8);
	for (size_t t = 0; t < totFrames; t++) {
		for (int j = 0; j < ini.totNbin; j++)
		{
			for (int i = 0; i < ini.ngrps; i++)
			{
				ofile << std::setw(15) << massdens[t][i][j] << "\t"
					<< std::setw(15) << numdensM[t][i][j] << "\t"
					<< std::setw(15) << chargedens[t][i][j] << "\t"
					<< std::setw(15) << numdensG[t][i][j] << "\t";
			}
			ofile << "\n";
		}
	}
	ofile.close();
	printf("Analysis done.\n");
}

OnflyPostHandle::OnflyPostHandle(std::string fnm)
{
	fid.open(fnm, std::ios::binary);
	if (!fid)
	{
		printf("Cannot open file \"%s\", exit!\n", fnm.c_str());
		std::exit(-1);
	}
}

void OnflyPostHandle::readHead()
{
	int magicNumber = 0;
	fid.read((char*)&magicNumber, sizeof(int));
	//std::cout<<magicNumber;
	if (magicNumber == 20080513)
	{
		fid.read((char*)&headSize, sizeof(int));
		int commandLineSize = 0;
		fid.read((char*)&commandLineSize, sizeof(int));
		commandLine.resize(commandLineSize);
		fid.read((char*)commandLine.data(), commandLine.size());
		fid.read((char*)&ngrps, sizeof(int));
		fid.read((char*)&nRegion, sizeof(int));
		fid.read((char*)&userint1, sizeof(int));
		fid.read((char*)&userint2, sizeof(int));
		fid.read((char*)&dt, sizeof(double));
		lowPos.resize(3 * nRegion);
		upPos.resize(3 * nRegion);
		nbin.resize(3 * nRegion);
		dbin.resize(3 * nRegion);
		fid.read((char*)lowPos.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)upPos.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)dbin.data(), sizeof(real) * 3 * nRegion);
		fid.read((char*)nbin.data(), sizeof(int) * 3 * nRegion);
		fid.read((char*)Lbox, sizeof(real) * 3);
	}
	else
	{
		std::cout << "Error format, cannot read head!" << std::endl;
		std::exit(-1);
	}
	totNbin = 0;
	//printf("why?");
	fflush(stdout);
	for (int i = 0; i < nRegion; i++)
		totNbin += nbin[i * 3] * nbin[i * 3 + 1] * nbin[i * 3 + 2];
	massdens.resize(ngrps);
	chargedens.resize(ngrps);
	numdensM.resize(ngrps);
	numdensG.resize(ngrps);
	for (int i = 0; i < ngrps; i++)
	{
		massdens[i].resize(totNbin);
		chargedens[i].resize(totNbin);
		numdensM[i].resize(totNbin);
		numdensG[i].resize(totNbin);
	}
}

bool OnflyPostHandle::readNextFrame()
{
	int magicNumber = 0;
	while (true)
	{
		fid.read((char*)&magicNumber, sizeof(int));
printf("get output?%d\n",magicNumber);
					fflush(stdout);
		if (magicNumber == 20201210)
		{
			for (int i = 0; i < ngrps; i++)
			{
				fid.read((char*)(massdens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(numdensM[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(chargedens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(numdensG[i].data()), sizeof(real) * totNbin);
			}
			return true;
		}
		else if (magicNumber == 20080513)
		{
			int ignoreSize = 0;
			fid.read((char*)&ignoreSize, sizeof(int));
			fid.seekg(ignoreSize - 2 * sizeof(int), std::ios::cur);
			magicNumber = 0;
			continue;
		}
		else
		{
			return false;
		}
	}
}
