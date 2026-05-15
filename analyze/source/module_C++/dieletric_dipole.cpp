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

using std::string, std::vector;
typedef	vector<vector<float>>  itp_matrix;
typedef	vector<vector<vector<float>>>  itp_matrix_3D;

/****************************************************************************/
/* This program calculates the Electric field intensity  in selected index group, and it is dynamic */
/* zengliang eidt 2019.1.9                              */
/*refere enetgy-peng(fengjie peng January 2018)                              */
/****************************************************************************/

static int nz = 1000;
static int num2ave = 1;
static real zmin = 20.0;
static real zmax = 30.0;
static int dipoleBin = 90;
const char* title = "Calculate the energy inside the molecule";
const string xaxis = "vertical-for time/molecule; level-for zbin", yaxis = "E(kJ/mol)";

double pbc(double dx,double box);
void dw_boundary(real* Lbox);

static void do_vcomponent(const char *traj_file, t_topology *top, const char *ndx_file,
					string filename, gmx_bool bMOL, const gmx_output_env_t* oenv, const char* tpr_file)
{
	std::cout << " \n get start \n ";
	t_trxframe  fr;
	t_trxstatus *status;
	char**	grpname;
	int		flags = TRX_READ_X;   /* read only position */
	int		step, grpNum = 1;
	int*	tot_atom_sel, **index;
	real	dbinz;
	real	Lbox[3] = { 0 };	

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
	
	// Open and initialize the xvg file
	FILE * outfile = xvgropen(filename.c_str(), title, xaxis, yaxis, oenv);
	step = 0;
	vector<double> charge_pro(nz, 0.0);
	vector<double> dipole_z(nz,0.0);

	do
	{
		charge_pro.assign(nz, 0.0);
		for (int i=0; i<tot_atom_sel[0]; i++){
			float nowz = fr.x[index[0][i]][2];
			if(nowz <= zmin || nowz >= zmax) continue;
			float charge = top->atoms.atom[index[0][i]].q;
			int me = (nowz - zmin)/dbinz;
			if(me < 0 || me >= nz) continue;
			float dV = Lbox[0] * Lbox[1] * dbinz;
			charge_pro[me] += charge / dV;
		}

		double cumulative_charge = 0.0;
		for (int i = 0; i < nz; i++) {
			cumulative_charge += charge_pro[i]; 
			dipole_z[i] += -cumulative_charge * dbinz;
		}

		for(int i = 0; i < nz; i++){
			fprintf(outfile, "%lf	", dipole_z[i]);
			dipole_z[i] = 0.0;
		}
		fprintf(outfile, "\n");

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
	static gmx_bool bMOL = FALSE;
	t_pargs pa[] =
	{
		{
			"-mol", FALSE, etBOOL, {&bMOL},
			"count by molecules or atoms, default is treated as atoms(it is better set to molecules)"
		},
		{
			"-sl", FALSE, etINT, {&nz},
			"Divide the box in #nr slices in z direction,set 1000"
		},
		{
			"-zmin", FALSE, etREAL, {&zmin},
			"min z to calculate (nm),set to 0 when no limit"
		},
		{
			"-zmax", FALSE, etREAL, {&zmax},
			"max z to calculate (nm),set to 0 when no limit"
		},
	};

	t_filenm fnm[] =
	{
		{ efTRX, NULL, NULL,  ffREAD },
		{ efTPS, NULL, NULL,  ffREAD },
		{ efNDX, NULL, NULL,  ffOPTRD },
		{ efXVG, NULL, "energy", ffWRITE },
	};
	#define NFILE asize(fnm)
	if(!parse_common_args(&argc,argv,
	                  PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT | PCA_CAN_TIME,
	                  NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv))	exit(0);

	/* Get the file name of the input/output file */
	traj_file	= ftp2fn_null(efTRX,NFILE,fnm);
	tpr_file	= ftp2fn_null(efTPS,NFILE,fnm);
	ndx_file	= ftp2fn_null(efNDX,NFILE,fnm);
	filename	= ftp2fn_null(efXVG,NFILE,fnm);

	/* read topology file，top contains almost all the information from the Force field file and the gro file*/
	top = read_top(tpr_file, &ePBC);

	do_vcomponent(traj_file, top, ndx_file, filename, bMOL, oenv, tpr_file);
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