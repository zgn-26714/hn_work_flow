/* yeting
 * 2020.12.14
 * yeting2938@hust.edu.cn
 * for onfly 3D density post ananysis
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "./include/core_getopt.h"


using real = double;
const std::vector<std::string> matlabcode = { {R"(% function [ret, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(fnm)
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
% end)"},
	{},
	{}
};

typedef std::vector<std::vector<real>> realMatrix;
typedef std::vector<std::vector<double>> doubleMatrix;

struct OnflyPostHandle
{
public:
	OnflyPostHandle(std::string fnm);
	void readHead();
	int readNextFrame();

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
	realMatrix massdens;
	realMatrix numdensM;
	realMatrix chargedens;
	realMatrix numdensG;
	realMatrix calcbond;
	realMatrix calcangle;
	doubleMatrix calcLJ;
	std::vector<double> allLJ;
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
	std::vector<std::string> outputFnms;
	double beginTime = 0;
	double endTime = 0;
	double numFile = 1;
	int dynamic = 1;
	int output_mode = 0;
	std::string outputFnm = "./deal_data/onfly/onflyEnergy";

	itp::Getopt getopt(argc, argv, "Onfly PostAnalysis Code");
	getopt.getFixPos(inputFnm, 1, true, "input file name");
	getopt(output_mode, "-om", false, "output_mode (ps)");
	getopt(beginTime, "-b", false, "begin time (ps)");
	getopt(endTime, "-e", false, "end time (ps)");
	getopt(outputFnm, "-o", false, "output file name");
	getopt(numFile, "-n", false, "input file number");
	getopt(dynamic, "-d", false, "isDynamic?");
	getopt.finish();

	outputFnms.emplace_back(outputFnm + "_dens.dat");
	outputFnms.emplace_back(outputFnm + "_BA.dat");
	outputFnms.emplace_back(outputFnm + "_LJ.dat");

	std::cout << "\nMake sure you have this stucture:"
		<< "\n\t./case1/file"
		<< "\n\t./case2/file"
		<< "\n\t./case3/file"
		<< "\n\t\t...";
	printf("beginTime: %f\n", beginTime);
	printf("endTine: %f\n", endTime);
	printf("input file name: %s\n", inputFnm.c_str());
	printf("output file name: %s\t%s\t%s\n", outputFnms[0].c_str(), outputFnms[1].c_str(), outputFnms[2].c_str());
	printf("input file number: %f\n", numFile);
	printf("isdynamic: %d\n", dynamic);

	OnflyPostHandle ini(inputFnm);
	ini.readHead();

	double dt = ini.dt * ini.userint2;
	size_t beginStep = beginTime / dt;
	size_t endStep = endTime / dt;
	size_t totFrames = 0;
	size_t timeindex;
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
//通过识别二进制文件中一帧结束的标识符来统计帧数
	for (size_t i = beginStep; i <= endStep; i++)
	{
		int index = ini.readNextFrame();
		if (!index)
			break;
		else if (100 == index) {
			totFrames++;
	//		std::cout<<totFrames<<'\n';fflush(stdout);
		}
	}
	timeindex = totFrames;
	if (0 == dynamic) totFrames = 1;
	std::cout << totFrames << "\n";
	doubleMatrix	allLJ;
	std::vector<realMatrix> massdens;
	std::vector<realMatrix> numdensM;
	std::vector<realMatrix> chargedens;
	std::vector<realMatrix> numdensG;
	std::vector<realMatrix> calcbond;
	std::vector<realMatrix> calcangle;
	std::vector<doubleMatrix> calcLJ;
	massdens.resize(totFrames);
	chargedens.resize(totFrames);
	numdensM.resize(totFrames);
	numdensG.resize(totFrames);
	calcbond.resize(totFrames);
	calcangle.resize(totFrames);
	calcLJ.resize(totFrames);
	allLJ.resize(totFrames);
	for (int i = 0; i < totFrames; i++){
		massdens[i].resize(ini.ngrps);
		chargedens[i].resize(ini.ngrps);
		numdensM[i].resize(ini.ngrps);
		numdensG[i].resize(ini.ngrps);
		calcbond[i].resize(ini.ngrps);
		calcangle[i].resize(ini.ngrps);
		calcLJ[i].resize((ini.ngrps + 3) * (ini.ngrps + 2) / 2);
		allLJ[i].resize(ini.nbin[2], 0.0);
		for (int j = 0; j < ini.ngrps; j++){
			massdens[i][j].resize(ini.totNbin, 0);
			numdensM[i][j].resize(ini.totNbin, 0);
			chargedens[i][j].resize(ini.totNbin, 0);
			numdensG[i][j].resize(ini.totNbin, 0);
			calcbond[i][j].resize(ini.totNbin,0.0);
			calcangle[i][j].resize(ini.totNbin,0.0);
		}
		for (int j = 0; j < calcLJ[0].size(); j++) {
			calcLJ[i][j].resize(ini.nbin[2], 0.0);
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
		if (1 == dynamic)
		{
			for (size_t t = 0; t < timeindex; t++)
			{
				while(true){
					int index = hd.readNextFrame();
					if (!index)
						break;
					else if (1 == index) {
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
					}
					else if (2 == index) {
						for (int i = 0; i < hd.ngrps; i++)
						{
							for (int j = 0; j < hd.totNbin; j++)
							{
								calcbond[t][i][j] += hd.calcbond[i][j];
								calcangle[t][i][j] += hd.calcangle[i][j];
							}
						}
					}
					else if (3 == index) {
						for (int j = 0; j < hd.nbin[2]; j++) {
							allLJ[t][j] += hd.allLJ[j];
							for (int i = 0; i < calcLJ[0].size(); i++) {
								calcLJ[t][i][j] += hd.calcLJ[i][j];
							}
						}
					}
					else if(100 == index) break;
				}
			}
		}
		else if (0 == dynamic) {
			std::cout << "TotFrames: " << nowTime << '\n';
			for (size_t t = 0; t < timeindex; t++)
			{
				while(true){
					int index = hd.readNextFrame();
					if (!index)
						break;
					else if (1 == index) {
						for (int i = 0; i < hd.ngrps; i++)
						{
							for (int j = 0; j < hd.totNbin; j++)
							{
								massdens[0][i][j] += hd.massdens[i][j];
								numdensM[0][i][j] += hd.numdensM[i][j];
								chargedens[0][i][j] += hd.chargedens[i][j];
								numdensG[0][i][j] += hd.numdensG[i][j];
							}
						}
					}
					else if (2 == index) {
						for (int i = 0; i < hd.ngrps; i++)
						{
							for (int j = 0; j < hd.totNbin; j++)
							{
								calcbond[0][i][j] += hd.calcbond[i][j];
								calcangle[0][i][j] += hd.calcangle[i][j];
							}
						}
					}
					else if (3 == index) {
						for (int j = 0; j < hd.nbin[2]; j++) {
							allLJ[0][j] += hd.allLJ[j];
							for (int i = 0; i < calcLJ[0].size(); i++) {
								calcLJ[0][i][j] += hd.calcLJ[i][j];
							}
						}
					}
					else if(100 == index) break;
				}                
			}
		}
		inputFnm = get_inputFile(inputFnm);
	}
	size_t averageNum = numFile;
	if (0 == dynamic) averageNum *= timeindex;
	for (size_t t = 0; t < totFrames; t++)
	{
		for (int i = 0; i < ini.ngrps; i++)
		{
			for (int j = 0; j < ini.totNbin; j++)
			{
				massdens[t][i][j] /= averageNum;
				numdensM[t][i][j] /= averageNum;
				chargedens[t][i][j] /= averageNum;
				numdensG[t][i][j] /= averageNum;
				calcbond[t][i][j] /= averageNum;
				calcangle[t][i][j] /= averageNum;
			}
		}
		for (int j = 0; j < ini.nbin[2]; j++) {
			allLJ[t][j] /= averageNum;
			for (int i = 0; i < calcLJ[0].size(); i++) {
				calcLJ[t][i][j] /= averageNum;
			}
		}
	}

	std::vector<std::ofstream> ofile(3);
	for (int file_idnex = 0; file_idnex < outputFnms.size(); file_idnex++)
	{
		ofile[file_idnex].open(outputFnms[file_idnex]);
		ofile[file_idnex] << "% onflyCommandline = \"" << ini.commandLine << "\";\n";
		ofile[file_idnex] << "% Lbox = [";
		for (auto&& i : ini.Lbox)
			ofile[file_idnex] << i << " ";
		ofile[file_idnex] << "];\n";
		ofile[file_idnex] << "% lowPos = [";
		for (auto&& i : ini.lowPos)
			ofile[file_idnex] << i << " ";
		ofile[file_idnex] << "];\n";
		ofile[file_idnex] << "% upPos = [";
		for (auto&& i : ini.upPos)
			ofile[file_idnex] << i << " ";
		ofile[file_idnex] << "];\n";
		ofile[file_idnex] << "% fnm = \"" << outputFnm << "\";\n";
		ofile[file_idnex] << "% nbin = [";
		for (auto&& i : ini.nbin)
		{
			ofile[file_idnex] << i << " ";
		}
		ofile[file_idnex] << "];\n";
		ofile[file_idnex] << "% nRegion = " << ini.nRegion << ";\n";
		ofile[file_idnex] << "% totNbin = " << ini.totNbin << ";\n";
		ofile[file_idnex] << "% totframes = " << totFrames << ";\n";
		ofile[file_idnex] << "% [data, xyzRange, nbin, lowPos, upPos, Lbox] = loadOnflyData3D(fnm);\n";
		ofile[file_idnex] << "%\n";
		ofile[file_idnex] << matlabcode[file_idnex];
		ofile[file_idnex] << "\n";
		ofile[file_idnex].setf(std::ios::scientific);
		ofile[file_idnex].precision(8);
	}

	for (size_t t = 0; t < totFrames; t++) {
		for (int j = 0; j < ini.totNbin; j++)
		{
			for (int i = 0; i < ini.ngrps; i++)
			{
				ofile[0] << std::setw(15) << massdens[t][i][j] << "\t"
					<< std::setw(15) << numdensM[t][i][j] << "\t"
					<< std::setw(15) << chargedens[t][i][j] << "\t"
					<< std::setw(15) << numdensG[t][i][j] << "\t";
			}
			ofile[0] << "\n";
		}
	}
	ofile[0].close();

	if (output_mode != 1 && output_mode != 3) {
		for (size_t t = 0; t < totFrames; t++) {
			for (int j = 0; j < ini.totNbin; j++)
			{
				for (int i = 0; i < ini.ngrps; i++)
				{
					ofile[1] << std::setw(15) << calcbond[t][i][j] << "\t"
						<< std::setw(15) << calcangle[t][i][j] << "\t";
				}
				ofile[1] << "\n";
			}
		}
		ofile[1].close();
	}
	
	if (output_mode != 2 && output_mode != 3) {
		for (size_t t = 0; t < totFrames; t++) {
			for (int j = 0; j < ini.nbin[2]; j++)
			{
				ofile[2] << std::setw(15) << allLJ[t][j] << '\t';
				for (int i = 0; i < (ini.ngrps + 3) * (ini.ngrps + 2) / 2; i++)
				{
					ofile[2] << std::setw(15) << calcLJ[t][i][j] <<'\t';
				}
				ofile[2] << "\n";
			}
		}
		ofile[2].close();
	}



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
	calcbond.resize(ngrps);
	calcangle.resize(ngrps);
	calcLJ.resize((ngrps + 3) * (ngrps + 2) / 2);
	allLJ.resize(nbin[2]);
	for (int i = 0; i < ngrps; i++)
	{
		massdens[i].resize(totNbin);
		chargedens[i].resize(totNbin);
		numdensM[i].resize(totNbin);
		numdensG[i].resize(totNbin);
		calcbond[i].resize(totNbin);
		calcangle[i].resize(totNbin);
	}
	for (int i = 0; i < (ngrps + 2) * (ngrps + 3) / 2; i++) {
		calcLJ[i].resize(nbin[2]);
	}
}

int OnflyPostHandle::readNextFrame()
{
	int magicNumber = 0;
	while (true)
	{
		fid.read((char*)&magicNumber, sizeof(int));
//printf("get output?%d\n",magicNumber);
//					fflush(stdout);
		if (magicNumber == 20201210)
		{
			for (int i = 0; i < ngrps; i++)
			{
				fid.read((char*)(massdens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(numdensM[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(chargedens[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(numdensG[i].data()), sizeof(real) * totNbin);
			}
			return 1;
		}
		else if (magicNumber == 20250416) {
			for (int i = 0; i < ngrps; i++)
			{
				fid.read((char*)(calcbond[i].data()), sizeof(real) * totNbin);
				fid.read((char*)(calcangle[i].data()), sizeof(real) * totNbin);
			}
			return 2;
		}
		else if (magicNumber == 20250417) {
			for (int i = 0; i < (ngrps + 2) * (ngrps + 3) / 2; i++) {
				fid.read((char*)(calcLJ[i].data()), sizeof(double) * nbin[2]);
			}
			fid.read((char*)(allLJ.data()), sizeof(double) * nbin[2]);
			return 3;
		}
		else if (magicNumber == 20080513)
		{
			int ignoreSize = 0;
			fid.read((char*)&ignoreSize, sizeof(int));
			fid.seekg(ignoreSize - 2 * sizeof(int), std::ios::cur);
			magicNumber = 0;
			continue;
		}
                else if(magicNumber == 20250429)
                {//end flag
                        return 100;
                }
		else
		{
			return 0;
		}
	}
}
