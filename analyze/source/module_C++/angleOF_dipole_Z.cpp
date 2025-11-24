#include <gromacs/fileio/confio.h>
#include <gromacs/fileio/xvgr.h>
#include <gromacs/utility/arraysize.h>
#include <gromacs/commandline/cmdlineinit.h>
#include <gromacs/trajectory/trajectoryframe.h>
#include <gromacs/utility/smalloc.h>
#include <gromacs/topology/topology.h>
#include <gromacs/topology/index.h>
#include <gromacs/utility/fatalerror.h>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/fileio/tpxio.h>
#include <gromacs/fileio/trxio.h>
#include <gromacs/commandline/pargs.h>

#ifdef HAVE_CONFIG_H
#endif
#include <numeric>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <cfloat>
#include "./include/matrix.h"
#include <cmath>



using std::string, std::vector;
using namespace itp;
const int OPEN		= 1;
const int CLOSE		= 0;
/****************************************************************************/
/* This program calculates the Electric field intensity  in selected index group, and it is dynamic */
/* zengliang eidt 2019.1.9                              */
/*refere enetgy-peng(fengjie peng January 2018)                              */
/****************************************************************************/

static int nz = 1000;
static int  rmind = 0;  /* this is for mismatch in index when using xtc file while xtc did not storage all groups*/
static int num2ave = 1;
static real zmin = 20.0;
static real zmax = 30.0;
static int  dipoleBin = 90;
static int out_sel=0;
const char* title = "Calculate the angle of dipole array and z array";
const string xaxis = "vertical-for time", yaxis = "dipole bin";

double pbc(double dx,double box);
void loadMessageAtom(Matrix_3d<real>& data, t_trxframe* fr, int* index, 
					 int N_mol, int Natom_mol, int message, rvec* f_global = NULL);
void loadMessageMol(Matrix<real>& posc, int com, t_trxframe* fr, int* index, 
					int N_mol, int Natom_mol, int message, vector<real>& mass,
					vector<real>& charge, rvec* f_global = NULL);
void dw_boundary(real* Lbox);
double array_angle(vector<double>& array1,vector<double>& array2);

static void do_vcomponent(	const char *traj_file, t_topology *top, const char *ndx_file,
					string filename, gmx_bool bMOL, const gmx_output_env_t* oenv, const char* tpr_file)
{
	std::cout << " \n get start \n ";
	t_trxframe  fr;
	t_trxstatus *status;
	char**	grpname;
	int		flags = TRX_READ_X;   /* read only position */
	int		step, nowz, grpNum = 1;
	int*	tot_atom_sel, **index;
	real	dbinz;
	real	Lbox[3] = { 0 };
	vector<int>		Natom_mol(grpNum, 0.0), N_mol(grpNum, 0.0);
	Matrix<real>	mass(grpNum, 0), charge(grpNum, 0);

	
	/* The first frame is read and the input boundary conditions are processed according to the box size */
	read_first_frame(oenv, &status, traj_file, &fr, flags);
	for (int i = 0; i < 3; i++) Lbox[i] = fr.box[i][i];
	dw_boundary(Lbox);
	dbinz = (zmax - zmin) / nz;
	std::cout << "\n now silce is nz=" << nz << '\n';
	/* Allocate memory and get index files */
	snew(index, grpNum);
	snew(tot_atom_sel, grpNum);
	snew(grpname, grpNum);
	get_index(&top->atoms, ndx_file, grpNum, tot_atom_sel, index, grpname);

	/* Get the number/mass/charge of atoms of each molecule */
	for (int i = 0; i < grpNum; i++) {
		int molindex = top->atoms.atom[index[i][0]].resind;
		string tmp = string(top->atoms.resinfo[molindex].name[0]);//这样写应该输入分子索引
		int count = 0;
		while (true) {
			if (top->atoms.atom[index[i][Natom_mol[i]]].resind == top->atoms.atom[index[i][0]].resind) {
				mass[i].emplace_back(top->atoms.atom[index[i][Natom_mol[i]]].m);
				charge[i].emplace_back(top->atoms.atom[index[i][Natom_mol[i]]].q);
				Natom_mol[i]++;
			}
			else	break;
		}
		N_mol[i] = tot_atom_sel[i] / Natom_mol[i];
	}
	
	// Open and initialize the xvg file
	FILE * outfile = xvgropen(filename.c_str(), title, xaxis, yaxis, oenv);
	step = 0;
	//auto start = std::chrono::high_resolution_clock::now();
    Matrix<real>    angle(grpNum, nz, 0);
	Matrix_3d<real>	angle_pro(grpNum, nz, dipoleBin, 0);
	Matrix<int>		count(grpNum, nz, 0);
	vector<double>	refArray = {0.0, 0.0, 1.0}; // reference array for angle calculation
	std::ofstream debug_file("./result/debug_angleOFdipoleZ.log");
	do
	{
		for (int i = 0; i < grpNum; i++) {
			Matrix_3d<real>		pos(N_mol[i], Natom_mol[i], 3, 0);
			Matrix<real>		posc(N_mol[i], 3, 0);
			loadMessageAtom(pos, &fr, index[i], N_mol[i], Natom_mol[i], 0);
			loadMessageMol(posc, 0, &fr, index[i], N_mol[i], Natom_mol[i], 0, mass[i], charge[i]);
			//#pragma omp parallel for num_threads(8)
			for (int j = 0; j != N_mol[i]; j++) {

				std::vector<double> ang_dipole(3, 0.0f);
				/* transform the molecule index to position index */
				if (zmin <= posc[j][2] && posc[j][2] < zmax) {
					if(posc[j][2] == zmax)	std::cout<<zmax<<std::endl;
					nowz = int((posc[j][2] - zmin) / dbinz);
					if (nowz < 0 || nowz >= nz) {
						debug_file << "[WARN] frame " << step
								<< " group " << i << " molecule " << j
								<< " nowz out of range: nowz=" << nowz
								<< " posz=" << posc[j][2]
								<< " range=[0," << nz-1 << "]" << std::endl;
						continue;
					}
					//getRefArray(pos[j], refArray);
					for (int num = 0; num < Natom_mol[i]; num++) {
						for(int k = 0; k < 3; k++){
							ang_dipole[k] += (pos[j][num][k] - posc[j][k]) * charge[i][num];
						}
					}
					double tmp_ang = array_angle(ang_dipole, refArray);
					if (!std::isfinite(tmp_ang)) {
						debug_file << "[WARN] frame " << step
								<< " group " << i << " molecule " << j
								<< " invalid tmp_ang = " << tmp_ang << std::endl;
						continue;
					}
                    angle[i][nowz] += tmp_ang;
					int bin = int(tmp_ang / (180.0 / dipoleBin)); // 180.0 is the max angle
					if (bin < 0) bin = 0;
                    if (bin >= dipoleBin) {
						debug_file << "[WARN] frame " << step
								<< " group " << i << " molecule " << j
								<< " bin out of range: bin=" << bin
								<< " tmp_ang=" << tmp_ang
								<< " clamp to " << dipoleBin-1 << std::endl;
						bin = dipoleBin - 1;
					}
					// if(bin == 180)	std::cout<<180<<std::endl;
					angle_pro[i][nowz][bin] ++;
					count[i][nowz] ++;
				}
				
			}
		}
		if (step%num2ave == (num2ave-1)){
			//#pragma omp parallel for num_threads(32)
			for (int i = 0; i < grpNum; i++) {
				for (int j = 0; j < nz; j++) {
                    if(out_sel){
                        if(count[i][j] != 0)
                            angle[i][j] /=  count[i][j];
                        fprintf(outfile, "%10g	", angle[i][j]);
                    }
                    else{
                        angle[i][j] = 0;
                        for (int k = 0; k < dipoleBin; k++) {
                            if(count[i][j] != 0)
                                angle_pro[i][j][k] /= count[i][j];
                            fprintf(outfile, "%10g	", angle_pro[i][j][k]);
                            angle_pro[i][j][k] = 0;
                            
					    }
                    }
					count[i][j] = 0;
                    angle[i][j] = 0;
                    std::fill(angle_pro[i][j].begin(),angle_pro[i][j].end(), 0);
					fprintf(outfile, "\n");
				}
			}
		}
		step++;
	}while(read_next_frame(oenv,status,&fr));
	fclose(outfile);
}

int main_func(int argc,char *argv[])
{
	t_topology*		top;
	int				ePBC;
	const char*		traj_file, * tpr_file, * ndx_file, * filename;
	gmx_output_env_t* oenv;

	const char *desc[] ={
		"This program calculates the h-component-v, include dipole and quadrupole.",
	};
	static gmx_bool bMOL       = FALSE;
	t_pargs pa[] =							//The customization parameters of the program
	{
		{
			"-mol",    FALSE, etBOOL, {&bMOL},
			"count by molecules or atoms, default is treated as atoms(it is better set to molecules)"
		},
		{
			"-rm",FALSE, etINT, {&rmind},
			"number of atoms didn't storage in xtc file before selected group"
		},
		{
			"-n2ave",FALSE, etINT, {&num2ave},
			"average of time,set 1"
		},
		{
			"-sl",FALSE, etINT, {&nz},
			"Divide the box in #nr slices in z direction,set 1000"
		},
		{
			"-zmin",FALSE, etREAL, {&zmin},
			"min z to calculate (nm),set to 0 when no limit"
		},
		{
			"-zmax",FALSE, etREAL, {&zmax},
			"max z to calculate (nm),set to 0 when no limit"
		},
		{
			"-dBin",FALSE, etINT, {&dipoleBin},
			"max z to calculate (nm),set to 0 when no limit"
		},
        {
			"-out",FALSE, etINT, {&out_sel},
			"max z to calculate (nm),set to 0 when no limit"
		}
	};

	t_filenm fnm[] =						//Necessary input/output files
	{
		{ efTRX, NULL, NULL,  ffREAD },
		{ efTPS, NULL, NULL,  ffREAD },
		{ efNDX, NULL, NULL,  ffOPTRD },
		{ efXVG, NULL, "angleOFdz", ffWRITE },
	};
	#define NFILE asize(fnm)				//Define the number of files as macros
	if(!parse_common_args(&argc,argv,		//Parsing command line
	                  PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT | PCA_CAN_TIME,
	                  NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv))	exit(0);

	/* Get the file name of the input/output file */
	traj_file	= ftp2fn_null(efTRX,NFILE,fnm);
	tpr_file	= ftp2fn_null(efTPS,NFILE,fnm);
	ndx_file	= ftp2fn_null(efNDX,NFILE,fnm);
	filename	= ftp2fn_null(efXVG,NFILE,fnm);

	/* read topology file，top contains almost all the information from the Force field file and the gro file*/
	top = read_top(tpr_file, &ePBC);

	/* Reading the bond and angle parameters is not necessary */
	do_vcomponent(traj_file, top, ndx_file, filename, bMOL, oenv, tpr_file);
	/*  the function end  */
	return 0;
}
double pbc(double dx,double box)
{
	while (dx > box/2.0)
		dx -=  box;
	while (dx < -box/2.0)
		dx +=  box;
	return dx;
}

void loadMessageMol(Matrix<real>& posc, int com, t_trxframe* fr, int* index,
					int N_mol, int Natom_mol, int message, vector<real>& mass,
					vector<real>& charge, rvec* f_global)
{
	int    i, j, k, molN;
	int    n = 0;
	real tmp, tmpCenter, halfLbox[3] = { 0 };
	vector<real> atomCenter;
	real         centerSum = 0;

	if (com == 0)
	{
		atomCenter = mass;
	}
	else if (com == 1)
	{
		atomCenter.resize(Natom_mol);
		std::fill_n(atomCenter.begin(), atomCenter.size(), 1);
	}
	else if (com == 2)
	{
		atomCenter = charge;
	}
	centerSum = std::accumulate(atomCenter.begin(), atomCenter.end(), real(0.0));

	for (k = 0; k < 3; k++)
		halfLbox[k] = fr->box[k][k] / 2;

	for (i = 0; i != N_mol; i++)
	{
		molN = index[n];
		for (k = 0; k != 3; k++)
		{
			tmpCenter = 0.0; // initiate for each molecule
			for (j = 0; j != Natom_mol; j++)
			{
				if (message == 0) {
					tmp = fr->x[index[n + j]][k];
					// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
					while (tmp - fr->x[molN][k] > halfLbox[k])
						tmp -= fr->box[k][k];
					while (tmp - fr->x[molN][k] < -halfLbox[k])
						tmp += fr->box[k][k];
				}
				else if (message == 1) {
					tmp = fr->v[index[n + j]][k];
				}
				else if (message == 2) {
					tmp = f_global[index[n + j]][k];
				}
				tmpCenter += tmp * atomCenter[j];
			}
			posc[i][k] = tmpCenter / centerSum;
		}
		n += Natom_mol;
	}
}

void loadMessageAtom(Matrix_3d<real>& data,  t_trxframe* fr, int* index,
					 int N_mol,  int Natom_mol,int message, rvec* f_global)			 
{
	int    i, j, k;
	int    molN;
	int    n = 0;
	real tmp[3]{}, halfLbox[3]{};

	for (k = 0; k < 3; k++)
		halfLbox[k] = fr->box[k][k] / 2;

	for (i = 0; i != N_mol; i++)
	{
		molN = index[n];
		for (j = 0; j != Natom_mol; j++)
		{
			for (k = 0; k != 3; k++)
			{
				if (message == 0)
				{
					tmp[k] = fr->x[index[n + j]][k];
					while (tmp[k] - fr->x[molN][k] > halfLbox[k])
						tmp[k] -= fr->box[k][k];
					while (tmp[k] - fr->x[molN][k] < -halfLbox[k])
						tmp[k] += fr->box[k][k];
				}
				else if (message == 1)
					tmp[k] = fr->v[index[n + j]][k];
				else if (message == 2)
					tmp[k] = f_global[index[n + j]][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--fr.x[molN][k] as absolute location)
				data[i][j][k] = tmp[k];
			}
		}
		n += Natom_mol;
	}
}


double array_angle(vector<double>& array1, vector<double>& array2)
{
	double up, down;
	long double pi = acos(0.0) * 2;
	up = array1[0] * array2[0] + array1[1] * array2[1] + array1[2] * array2[2];
	down =	sqrt(array1[0] * array1[0] + array1[1] * array1[1] + array1[2] * array1[2]) *
			sqrt(array2[0] * array2[0] + array2[1] * array2[1] + array2[2] * array2[2]);

	return (acos(up/down) * (180.0 / pi));//rad to deg
	//return (acos(up / down));
}

void dw_boundary(real *Lbox)
{
	if (zmin <= 0) {
		std::cout << "\nchange zmin from " << zmin << " to 0\n";
		zmin = 0;
	}
	if (zmax <= 0 || zmax > Lbox[2]) {
		std::cout << "\nchange zmax from " << zmax << " to "<<Lbox[2]<<'\n';
		zmax = Lbox[2];
	}
	if (zmax <= zmin) {
		std::cout << "\n\terror!	zmax must bigger than zmin!\n";
		exit(-1);
	}
}

int main(int argc, char** argv)
{
	return gmx_run_cmain(argc, argv, &main_func);
}
