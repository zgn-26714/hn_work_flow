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


/****************************************************************************/
/* This program calculates the Electric field intensity  in selected index group, and it is dynamic */
/* zengliang eidt 2019.1.9                              */
/*refere enetgy-peng(fengjie peng January 2018)                              */
/****************************************************************************/
static real xmin = 0;
static real ymin = 0;
static real zmin = 0;
static real xmax = 0;
static real ymax = 0;
static real zmax = 0;
static rvec nslice= {1, 1, 1000};
static int  rmind = 0;  /* this is for mismatch in index when using xtc file while xtc did not storage all groups*/
static int num2ave=1;
static rvec min   = {0.0, 0.0, 19.0};
static rvec max   = {4.5, 4.5, 29.0};
static real kB    = 1.3806505;/* Boltzmann constant without 10^-23 */
static real NA    = 6.02214129;  /* Avogadro constant without 10^23 */
static real fcomb = 138.935485; /* kJ mol^−1 nm e^−2 factor in coulomb interaction */
double * CreateVector(int rows);
double ** CreateMatrix(int rows,int cols);
double *** CreateMatrix_3d(int m,int n, int p);
double pbc(double dx,double box);
//next 3 function not use
int load_positionCenter(rvec *x0,int nm,int napm,int nStart,double **pos,double Lbox[3],double *deloc, double totMNC);
int load_positionCenter_fr(t_trxframe fr,int nm,int napm,int nStart,double **pos,double Lbox[3],double *deloc, double totMNC);
int load_position(rvec *x0,int nm,int napm,int nStart,double ***posmolecule,double Lbox[3]);

int load_position_force_fr(t_trxframe fr,int na,int **index1,double **pos1,double **force1,double Lbox[3]);
static void corr_print(const char *fn,const char *title,const char *yaxis,double dbinx,double dbiny,double dbinz,int nx,int ny,int nz,
                       double ***en1,double ***en2,double ***en3,char *grpname[],const gmx_output_env_t* oenv);

void do_forcebin(const char *trx_file, const char *ndx_file, const char *Ez_filenm,FILE *Ex_file,FILE *Ey_file,gmx_bool bMOL,
                 t_topology *top,const gmx_output_env_t* oenv)
{
	t_trxframe  fr;
	rvec        *x0;      /* coordinates without pbc */
	t_trxstatus *status;
	real        t,kden;
	const char  *ylabel = NULL;
	//int         natoms;  /* nr. atoms in trj */
	int         i,j,k,w,nx,ny,nz,mex,mey,mez,step,nm,nstart,napm,bbox[3];
	int         ii,k1,k2,type,Natom;
	int         flags = TRX_READ_X | TRX_READ_F;   /* read only position */
	double      Lbox[3],dt,halfLbox[3],*qi,dx,dy,dz,sumq,nowq;
	int     **index;
	int         *tot_atom; /* the selected groups' sizes */
	char        **grpname;
	double      **posCen,**pos1,**force1,dbinx,dbiny,dbinz,***Ex,***Ey,***Ez,***count;
	FILE *Ez_file;//dynamic outfile of Ez.Ex and Ey are opended in main funcion.
	FILE *outputfile=NULL;//tmp file

	read_first_frame(oenv,&status,trx_file,&fr,flags);
	//natoms=read_first_x(oenv, &status, trx_file, &t, &x0, box);
	Lbox[0] = fr.box[XX][XX];
	Lbox[1] = fr.box[YY][YY];
	Lbox[2] = fr.box[ZZ][ZZ];

	fprintf(stdout," \n get start  \n ");
	snew(index,1);
	snew(tot_atom,1);
	snew(grpname,1);
	get_index(&top->atoms,ndx_file,1,tot_atom,index,grpname);

	bbox[0] = 0;
	bbox[1] = 0;
	bbox[2] = 0;
	if(xmax == 0)
	{
		xmax = Lbox[0];
		if(xmin == 0)
		{
			bbox[0] = 1;
		}
	}
	if(ymax == 0)
	{
		ymax = Lbox[1];
		if(ymin == 0)
		{
			bbox[1] = 1;
		}
	}
	if(zmax == 0)
	{
		zmax = Lbox[2];
		if(zmin == 0)
		{
			bbox[2] = 1;
		}
	}
	nx   = (int)(nslice[0]+0.5);
	ny   = (int)(nslice[1]+0.5);
	nz   = (int)(nslice[2]+0.5);
	fprintf(stdout," \n now silce is \nx=%d \ny=%d\nz=%d\n ",nx,ny,nz);
	/*
	napm   = top->mols.index[1+top->atoms.atom[nstart].resind] - top->mols.index[top->atoms.atom[nstart].resind];
	if(!bMOL){napm = 1;}
	if (napm!=1) fprintf(stdout,"\n warning \n \n warning \n \n warning \n \n \
	warning \n \n warning \n the calculate is by center of molecular not atom \n");
	nm     = tot_atom[0]/napm;
	*/
	dbinx  = (xmax - xmin)/nx;
	dbiny  = (ymax - ymin)/ny;
	dbinz  = (zmax - zmin)/nz;
	pos1   = CreateMatrix(tot_atom[0],3);
	force1   = CreateMatrix(tot_atom[0],3);
	Ex = CreateMatrix_3d(nx,ny,nz);
	Ey = CreateMatrix_3d(nx,ny,nz);
	Ez = CreateMatrix_3d(nx,ny,nz);
	count  = CreateMatrix_3d(nx,ny,nz);
	for(i=0; i<nx; i++)
	{
		for(j=0; j<ny; j++)
		{
			for(k=0; k<nz; k++)
			{
				Ex[i][j][k] = 0.0;
				Ey[i][j][k] = 0.0;
				Ez[i][j][k] = 0.0;
				count[i][j][k]  = 0.0;
			}
		}
	}
	Natom = top->atoms.nr; // total atoms number of the system
	qi = CreateVector(Natom);
	for(i=0; i<Natom; i++)
	{
		qi[i] = top->atoms.atom[i].q;
	}// output first title and others
	Ez_file=xvgropen(Ez_filenm,"Electric filed(z-axis).",output_env_get_xvgr_tlabel(oenv),"E(kJ/mol/nm/e)",oenv);
	for(w=0;w<3 ;w++)
	{
		switch (w)
		{
			case 0:
			outputfile=Ez_file;break;
			case 1:
			outputfile=Ex_file;break;
			case 2:
			outputfile=Ey_file;break;
		}
		if (outputfile==NULL) continue;
		fprintf(outputfile,"# electric filed caculted frome group %5s \n",grpname[0]);
		fprintf(outputfile,"# calculate area: \n#  xmin-xmax  ymin-ymax  zmin-zmax:\n# number of frames to average is   %d\n",num2ave);
		fprintf(outputfile,"# calculate area: \n#  xmin-xmax  ymin-ymax  zmin-zmax:\n#  %.2f-%.2f, %.2f-%.2f, %.2f-%.2f,\n",xmin,xmax,ymin,ymax,zmin,zmax);
		fprintf(outputfile,"# the next line is location(or coordinate) x,y,z (nm) ,and the next next lines is Ez(kJ/mol/nm/e)  of different time  \n");
		for(i=0; i<nx; i++)
		for(j=0; j<ny; j++)
			for(k=0; k<nz; k++)
				fprintf(outputfile," %4g %4g %4g ",min[0]+(i+0.5)*dbinx,min[1]+(j+0.5)*dbiny,min[2]+(k+0.5)*dbinz);
				fprintf(outputfile,"\n");
				
	}

	step = 0;
	int step_m=0;
	do
	{
		load_position_force_fr(fr,tot_atom[0],index,pos1,force1,Lbox); // atom position and force for group1
		for(i=0; i<tot_atom[0]; i++)
		{
			mex = (int)((pos1[i][0]-xmin)/dbinx + 1) -1;
			mey = (int)((pos1[i][1]-ymin)/dbiny + 1) -1;
			mez = (int)((pos1[i][2]-zmin)/dbinz + 1) -1;
			if (mex>nx-1 || mex<0) continue;
			if (mey>ny-1 || mex<0) continue;
			if (mez>nz-1 || mex<0) continue;
			nowq=qi[index[0][i]];
			/*fprintf(stdout," \nstep is %d	%d\n",step,step_m);
			//if (mez>230 && mez<250){
			fprintf(stdout," \n number is %d	nowq %lf 	force %lf",index[0][i],nowq,force1[i][2]);
			}*/
			if (fabs(nowq)<1e-5) continue;
			count[mex][mey][mez]++;
			Ex[mex][mey][mez]+=force1[i][0]/nowq;
			Ey[mex][mey][mez]+=force1[i][1]/nowq;
			Ez[mex][mey][mez]+=force1[i][2]/nowq;


		}
		if(step%num2ave==(num2ave-1))//clean data and printf
			{
				for(i=0; i<nx; i++)
					for(j=0; j<ny; j++)
						for(k=0; k<nz; k++)
						{
							if(count[i][j][k]>0)
							{
								Ex[i][j][k]   /= count[i][j][k];
								Ey[i][j][k]   /= count[i][j][k];
								Ez[i][j][k]   /= count[i][j][k];
							}
							// output
							for(w=0;w<3;w++)
							{
								switch (w)
								{
									case 0:
									if (Ez_file!=NULL) fprintf(Ez_file,"%16g",Ez[i][j][k]);break; 
									case 1:
									if (Ex_file!=NULL) fprintf(Ex_file,"%16g",Ex[i][j][k]);break; 
									case 2:
									if (Ey_file!=NULL) fprintf(Ey_file,"%16g",Ex[i][j][k]);break;
								}									
							}
							//initial
							Ex[i][j][k]   = 0;
							Ey[i][j][k]   = 0;
							Ez[i][j][k]   = 0;
							count[i][j][k]= 0;
						}
						if(Ez_file!=NULL) fprintf(Ez_file,"\n");
						if(Ex_file!=NULL) fprintf(Ex_file,"\n");
						if(Ey_file!=NULL) fprintf(Ey_file,"\n");

			}
		step++;

		//fprintf(stdout," \n step=%d\n",step);

	}
	while(read_next_frame(oenv,status,&fr));
	fclose(Ez_file);
	if (Ex_file!=NULL)
	fclose(Ex_file);
	if (Ey_file!=NULL)
	fclose(Ey_file);

}

int main_func(int argc,char *argv[])
{
	const char *desc[] =
	{
		"This program calculates the force by divide \
	slices.",
	"force is electric force, input trr must rerun data which vdw forece is removed"
	};
	static gmx_bool bMOL       = FALSE;
	t_pargs pa[] =
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
			"number of atoms didn't storage in xtc file before selected group"
		},
		{
			"-sl",FALSE, etRVEC, {nslice},
			"Divide the box in #nr slices in x,y,z direction"
		},
		{
			"-min",FALSE, etRVEC, {min},
			"min (x,y,z) to calculate (nm),set to 0 when no limit"
		},
		{
			"-max",FALSE, etRVEC, {max},
			"max (x,y,z) to calculate (nm),set to 0 when no limit"
		}
	};

	t_filenm fnm[] =
	{
		{ efTRX, NULL, NULL,  ffREAD },
		{ efTPS, NULL, NULL,  ffREAD },
		{ efNDX, NULL, NULL,  ffOPTRD },
		{ efXVG, NULL, "Efield_z", ffWRITE },
		{ efXVG, "-ox", "Efield_x", ffOPTWR },
		{ efXVG, "-oy", "Efield_y", ffOPTWR },
	};
#define NFILE asize(fnm)

	t_topology  *top;
	int         ePBC;
	const char  *trx_file, *tps_file, *ndx_file, *Ez_filenm;
	FILE *Ex_file=NULL;
	FILE *Ey_file=NULL;
	gmx_bool bEx,bEy;
	//int         *gnx; /* the selected groups' sizes */
	//char        **grpname;
	gmx_output_env_t* oenv;


	if(!parse_common_args(&argc,argv,
	                  PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT | PCA_CAN_TIME,
	                  NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv))
{
	exit(0);
}
	trx_file = ftp2fn_null(efTRX,NFILE,fnm);
	tps_file = ftp2fn_null(efTPS,NFILE,fnm);
	ndx_file = ftp2fn_null(efNDX,NFILE,fnm);
	Ez_filenm = ftp2fn_null(efXVG,NFILE,fnm);
	
	bEx=opt2bSet("-ox", NFILE, fnm);
	bEy=opt2bSet("-oy", NFILE, fnm);
	if (bEx) Ex_file=xvgropen(opt2fn("-ox", NFILE, fnm), "Electric filed(x-axis).",
                              output_env_get_xvgr_tlabel(oenv), "E(kJ/mol/nm/e)", oenv);
	if (bEy) Ey_file=xvgropen(opt2fn("-oy", NFILE, fnm), "Electric filed(y-axis).",
                              output_env_get_xvgr_tlabel(oenv), "E(kJ/mol/nm/e)", oenv);
	top = read_top(tps_file,&ePBC);     /* read topology file */
// get_index(&top->atoms,ndx_file,1,gnx,index,grpname);
	/* the function you are going to use */
	xmin = min[0];
	ymin = min[1];
	zmin = min[2];
	xmax = max[0];
	ymax = max[1];
	zmax = max[2];
	do_forcebin(trx_file,ndx_file,Ez_filenm,Ex_file,Ey_file,bMOL,
	            top,oenv);
	/*  the function end  */

	return 0;
}

double * CreateVector(int rows)
{
	double  *m;

	m = (double *) calloc((unsigned int) rows,sizeof(double ));

	return m;
}
double ** CreateMatrix(int rows,int cols)
{
	int  i;
	double  **m;

	m = (double **) calloc((unsigned int) rows,sizeof(double *));
	for (i=0; i < rows; i++)
	{
		m[i] =(double *)  calloc((unsigned int) cols,sizeof(double ));
	}

	return m;
}
double *** CreateMatrix_3d(int m,int n, int p)
{
	int i;
	double ***m2;
	m2 = (double ***) calloc((unsigned int) m, sizeof(double **));
	for(i=0; i<m; i++)
		m2[i] = CreateMatrix(n,p);
	return m2;
}
double pbc(double dx,double box)
{
	while (dx > box/2.0)
		dx -=  box;
	while (dx < -box/2.0)
		dx +=  box;

	return dx;
}
int load_positionCenter_fr(t_trxframe fr ,int nm,int napm,int nStart,double **pos,double Lbox[3],double *deloc, double totMNC)
{
	//load_positionCenter_fr is edited by zengliang(2019.1.9), compared load_positionCenter input is fr not
	// x0,becasuse read_next_frame,not read_next_x
	int i,j,k,molN,m;
	double center[3], tmpPos,tmpCenter,halfLbox[3];
	for(k=0; k<3; k++)
		halfLbox[k] = Lbox[k]/2;

	for(i=0; i<nm; i++)
	{
		molN = nStart+i*napm;

		for(k=0; k<3; k++)
		{
			tmpCenter = 0.0;  // initiate for each molecule
			for(j=0; j<napm; j++)
			{
				tmpPos =fr.x[molN+j][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)
				while(tmpPos-fr.x[molN][k]> halfLbox[k])
					tmpPos -= Lbox[k];
				while(tmpPos-fr.x[molN][k]<-halfLbox[k])
					tmpPos += Lbox[k];
				tmpCenter += tmpPos*deloc[j];
			}
			pos[i][k] = tmpCenter/totMNC;
			// move molecules into box.
			while(pos[i][k]>Lbox[k])
				pos[i][k] -= Lbox[k];
			while(pos[i][k]<0)
				pos[i][k] += Lbox[k];
		}
	}
	return 0;
}
//int load_position_force_fr(t_trxframe fr,int na,int nStart,double pos1,double force1,double Lbox[3])
int load_position_force_fr(t_trxframe fr,int na,int **index1,double **pos1,double **force1,double Lbox[3])
{
	//load_positionCenter_fr is edited by zengliang(2019.1.9),load position and force atoms
	// x0,becasuse read_next_frame,not read_next_x
	int i,j,k,molN,m;
	double tmpPos;
	for(i=0; i<na; i++)
	{
		for(k=0; k<3; k++)
		{
			tmpPos =fr.x[index1[0][i]][k];
			while(tmpPos> Lbox[k])
				tmpPos -= Lbox[k];
			while(tmpPos< 0)
				tmpPos += Lbox[k];
			pos1[i][k] = tmpPos;
			force1[i][k]=fr.f[index1[0][i]][k];
		}
	}


	return 0;
}
int load_positionCenter(rvec *x0,int nm,int napm,int nStart,double **pos,double Lbox[3],double *deloc, double totMNC)
{
	int i,j,k,molN,m;
	double center[3], tmpPos,tmpCenter,halfLbox[3];
	for(k=0; k<3; k++)
		halfLbox[k] = Lbox[k]/2;

	for(i=0; i<nm; i++)
	{
		molN = nStart+i*napm;

		for(k=0; k<3; k++)
		{
			tmpCenter = 0.0;  // initiate for each molecule
			for(j=0; j<napm; j++)
			{
				tmpPos = x0[molN+j][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)
				while(tmpPos-x0[molN][k]> halfLbox[k])
					tmpPos -= Lbox[k];
				while(tmpPos-x0[molN][k]<-halfLbox[k])
					tmpPos += Lbox[k];
				tmpCenter += tmpPos*deloc[j];
			}
			pos[i][k] = tmpCenter/totMNC;
			// move molecules into box.
			while(pos[i][k]>Lbox[k])
				pos[i][k] -= Lbox[k];
			while(pos[i][k]<0)
				pos[i][k] += Lbox[k];
		}
	}

	return 0;
}
int load_position(rvec *x0,int nm,int napm,int nStart,double ***posmolecule,double Lbox[3])
{
	int i,j,k,molN;
	double center[3], tmpPos[3],halfLbox[3];
	for(k=0; k<3; k++)
		halfLbox[k] = Lbox[k]/2;

	for(i=0; i<nm; i++)
	{
		molN = nStart+i*napm;
		for(j=0; j<napm; j++)
			for(k=0; k<3; k++)
			{
				tmpPos[k] = x0[molN+j][k];
				// first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)
				while(tmpPos[k]-x0[molN][k]> halfLbox[k])
					tmpPos[k] -= Lbox[k];
				while(tmpPos[k]-x0[molN][k]<-halfLbox[k])
					tmpPos[k] += Lbox[k];

				posmolecule[i][j][k] = tmpPos[k];
			}
	}

	return 0;
}
static void corr_print(const char *fn,const char *title,const char *yaxis,double dbinx,double dbiny,double dbinz,int nx,int ny,int nz,
                       double ***en1,double ***en2,double ***en3,char *grpname[],const gmx_output_env_t* oenv)
{
	FILE *out;
	int  i,j,k;

	out=xvgropen(fn,title,output_env_get_xvgr_tlabel(oenv),yaxis,oenv);
	fprintf(out,"# electric filed caculted frome group %5s \n",grpname[0]);
	fprintf(out,"# calculate area: \n#  xmin-xmax  ymin-ymax  zmin-zmax:\n#  %.2f-%.2f, %.2f-%.2f, %.2f-%.2f,\n",xmin,xmax,ymin,ymax,zmin,zmax);
	fprintf(out,"# location x,y,z (nm),      Ex(kJ/mol/nm/e),   Ey(kJ/mol/nm/e),  Ez(kJ/mol/nm/e)    \n");
	// if(dens_opt[0][0] == 'n'){
	// fprintf(out,"# location x,y,z (nm), number density (#/nm^3)\n");
	// }
	// if(dens_opt[0][0] == 'm'){
	// fprintf(out,"# location x,y,z (nm), mass density (kg/m^3)\n");
	// }
	// if(dens_opt[0][0] == 'c'){
	// fprintf(out,"# location x,y,z (nm), charge density (e/nm^3)\n");
	// }

	for(i=0; i<nx; i++)
	{
		for(j=0; j<ny; j++)
		{
			for(k=0; k<nz; k++)
			{
				fprintf(out,"  %10g    ",min[0]+(i+0.5)*dbinx);
				fprintf(out,"%10g    ",min[1]+(j+0.5)*dbiny);
				fprintf(out,"%10g    ",min[2]+(k+0.5)*dbinz);
				fprintf(out,"%10g",en1[i][j][k]);
				fprintf(out,"      %10g",en2[i][j][k]);
				fprintf(out,"      %10g",en3[i][j][k]);
				fprintf(out,"\n");
			}
		}
	}
	fclose(out);
}

int main(int argc, char** argv)
{
	return gmx_run_cmain(argc, argv, &main_func);
}//num:argc ; string: **argv

