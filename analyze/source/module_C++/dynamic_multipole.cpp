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
#include <iostream>
#ifdef HAVE_CONFIG_H
#endif
using std::vector;

/****************************************************************************/
/* This program calculates the Electric field intensity  in selected index group, and it is dynamic */
/* zengliang eidt 2019.1.9                              */
/*refere enetgy-peng(fengjie peng January 2018)                              */
/****************************************************************************/
static real zmin = 0;
static real zmax = 0;
static int nz = 1000;
static int  rmind = 0;  /* this is for mismatch in index when using xtc file while xtc did not storage all groups*/
static int num2ave = 1;
static int ngrps = 3;
static double epsilon = 0.0000;//
double * CreateVector(int rows);
double ** CreateMatrix(int rows,int cols);
double *** CreateMatrix_3d(int m,int n, int p);
double pbc(double dx,double box);


int load_positionCenter_mass(t_trxframe fr ,int nm,int atoms,int *index,double **posc,double* mass,double Lbox[]);
int load_position(rvec *x0,int nm,int napm,int nStart,double ***posmolecule,double Lbox[]);
void hngetMolInfo(t_topology *top, int **index, int nm, int *atoms, double **mass, double **charge);
void getMolMultipole_z(t_trxframe fr, int atoms, int nm, int *index, double **posc, double Lbox[], double* divC, double *charge,
                                        double *dipole, double *quadrupole);
static void corr_print(const char *fn,const char *title,const char *yaxis,double dbinx,double dbiny,double dbinz,int nx,int ny,
                       double ***en1,double ***en2,double ***en3,char *grpname[],const gmx_output_env_t* oenv);

void do_multipole(const char *trx_file, const char *ndx_file, const char *file,gmx_bool bMOL,
                 t_topology *top,const gmx_output_env_t* oenv)
{
        t_trxframe  fr;
        rvec        *x0;      /* coordinates without pbc */
        t_trxstatus *status;
        real        t,kden;
        const char  *ylabel = NULL;
        //int         natoms;  /* nr. atoms in trj */
        int         nowz,step,bbox[3];
        int         ii,k1,k2,type,Natom;
        int         flags = TRX_READ_X;   /* read only position */
        double      Lbox[3],dbinz,halfLbox[3];
        int        **index;
        int         *tot_atom; /* the selected groups' sizes */
        char        **grpname;
        double      ***posc, **charge, **dipole, **quadrupole;
        double** divC, ** mass;
        double** Charge, ** Dipole_E, ** Quadru_E;
        FILE *outfile;
        read_first_frame(oenv,&status,trx_file,&fr,flags);
        //natoms=read_first_x(oenv, &status, trx_file, &t, &x0, box);
        Lbox[0] = fr.box[0][0];
        Lbox[1] = fr.box[1][1];
        Lbox[2] = fr.box[2][2];
        dbinz=(zmax-zmin)/nz;
        fprintf(stdout," \n get start \n ");


        posc = new double**[ngrps];
        charge = new double* [ngrps];
        dipole = new double* [ngrps];
        quadrupole = new double* [ngrps];
        divC = new double* [ngrps];
        mass = new double* [ngrps];

        snew(index, ngrps);
        snew(tot_atom, ngrps);
        snew(grpname, ngrps);//杩欓噷浼氬奖鍝嶅懡浠よ閲岃閫夊嚑涓垎缁?
        get_index(&top->atoms,ndx_file, ngrps,tot_atom,index,grpname);//selected group index
        bbox[2] = 0;
        if(zmax == 0)
        {
                zmax = Lbox[2];
                if(zmin == 0)
                {
                        bbox[2] = 1;
                }
        }
        fprintf(stdout," \n now silce is nz=%d\n ",nz);
        vector<int> atoms(ngrps,0), nm(ngrps,0);
        for (int grp = 0; grp != ngrps; grp++) {
                for (;;)
                {
                        if (top->atoms.atom[index[grp][atoms[grp]]].resind == top->atoms.atom[index[grp][0]].resind)
                                atoms[grp]++;//in this ,resind mean the molecule index
                        //top include something銆俿uch as atoms;atoms need some index ,every idnex include mass銆乧harge銆乤nd so on;
                        else
                                break;
                }
                nm[grp] = tot_atom[grp] / atoms[grp];
                posc[grp] = CreateMatrix(nm[grp], 3);
                charge[grp] = CreateVector(nm[grp]);
                dipole[grp] = CreateVector(nm[grp]);
                quadrupole[grp] = CreateVector(nm[grp]);
        }


        outfile=xvgropen(file,"Multipole expansion(z-axis).",output_env_get_xvgr_tlabel(oenv),"rho (e/nm^3)",oenv);
        Charge = CreateMatrix(ngrps, nz);
        Dipole_E = CreateMatrix(ngrps, nz);
        Quadru_E = CreateMatrix(ngrps, nz);

        for (int grp = 0; grp != ngrps; grp++)
        {
                mass[grp] = CreateVector(atoms[grp]);
                divC[grp] = CreateVector(atoms[grp]);
                for (int i = 0; i < atoms[grp]; i++) {
                        mass[grp][i] = top->atoms.atom[index[grp][i]].m;
                        divC[grp][i] = top->atoms.atom[index[grp][i]].q;
                }
        }

        std::cout << "get molecule information: \n";
        for (int grp = 0; grp < ngrps; grp++) {
                std::cout << "\tmolecule " << grp << "\n";
                for (int i = 0; i < atoms[grp]; i++) {
                        std::cout << "\t\tmass:\t" << mass[grp][i] << '\n';
                        std::cout << "\t\tcharge:\t" << divC[grp][i] << '\n';
                }
                std::cout << '\n';
        }


        step = 0;
        do
        {

                for (int grp = 0; grp != ngrps; grp++) {
                        load_positionCenter_mass(fr, nm[grp], atoms[grp], index[grp], posc[grp], mass[grp], Lbox);
                        getMolMultipole_z(fr, atoms[grp], nm[grp], index[grp], posc[grp], Lbox, divC[grp], charge[grp], dipole[grp], quadrupole[grp]);
                        vector<double> tmp_qua(nz, 0.0);
                        for (int i = 0; i < nm[grp]; i++) {
                                if (posc[grp][i][2] <= zmax && posc[grp][i][2] >= zmin) {
                                        nowz = int((posc[grp][i][2] - zmin) / dbinz);
                                        double tmp_dphi_P = dipole[grp][i];
                                        Charge[grp][nowz] += charge[grp][i];

                                        Dipole_E[grp][nowz] += tmp_dphi_P;
                                        tmp_qua[nowz] = quadrupole[grp][i];
                                }
                        }
                        for (int iz = 0; iz < nz - 1; iz++) {
                                Quadru_E[grp][iz] += (tmp_qua[iz + 1] - tmp_qua[iz]) / dbinz;
                        }
                }

                if (step%num2ave == (num2ave-1)){
                        for (int grp = 0; grp != ngrps; grp++) {
                                for (int i = 0; i < nz; i++) {
                                        fprintf(outfile, "%16g ", Charge[grp][i] / num2ave);
                                        Charge[grp][i] = 0.0;
                                }
                                for (int i = 0; i < nz; i++) {
                                        fprintf(outfile, "%16g ", Dipole_E[grp][i] / num2ave);
                                        Dipole_E[grp][i] = 0.0;
                                }
                                for (int i = 0; i < nz; i++) {
                                        fprintf(outfile, "%16g ", Quadru_E[grp][i] / num2ave);
                                        Quadru_E[grp][i] = 0.0;
                                }
                                fprintf(outfile, "\n");
                        }
                }
                step++;
        }
        while(read_next_frame(oenv,status,&fr));
        fclose(outfile);
}

int main_func(int argc,char *argv[])
{
        const char *desc[] =
        {
                "This program calculates the Multipole expansion, include diipole and quadrupole.",

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
                        "average of time"
                },
                {
                        "-sl",FALSE, etINT, {&nz},
                        "Divide the box in #nr slices in z direction"
                },
                {
                        "-min",FALSE, etREAL, {&zmin},
                        "min z to calculate (nm),set to 0 when no limit"
                },
                {
                        "-max",FALSE, etREAL, {&zmax},
                        "max z to calculate (nm),set to 0 when no limit"
                },
                {
                        "-ngrps",FALSE, etINT, {&ngrps},
                        "number of selection groups"
                }
        };

        t_filenm fnm[] =
        {
                { efTRX, NULL, NULL,  ffREAD },
                { efTPS, NULL, NULL,  ffREAD },
                { efNDX, NULL, NULL,  ffOPTRD },
                { efXVG, NULL, "Multipole", ffWRITE },
        };
#define NFILE asize(fnm)

        t_topology  *top;
        int         ePBC;
        const char  *trx_file, *tps_file, *ndx_file, *file;
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
        file = ftp2fn_null(efXVG,NFILE,fnm);
        top = read_top(tps_file,&ePBC);     /* read topology file */
// get_index(&top->atoms,ndx_file,1,gnx,index,grpname);
        /* the function you are going to use */
        do_multipole(trx_file,ndx_file,file,bMOL,top,oenv);
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
int load_positionCenter_O(t_trxframe fr ,int nm,int atoms,int *index,double **posc,double *mass,double Lbox[])
{
        int i,j,k,molN;
        double tmpPos;
        for(i=0; i<nm; i++)
        {
                molN = index[i*atoms];
                for(k=0; k<3; k++)
                {
                        tmpPos =fr.x[molN][k];
                        posc[i][k] = tmpPos;
                        // move molecules into box.
                }
        }
        return 0;
}

int load_positionCenter_mass(t_trxframe fr ,int nm,int atoms,int *index,double **posc,double* mass,double Lbox[])
{
        //load_positionCenter_fr is edited by zengliang(2019.1.9), compared load_positionCenter input is fr not^M
        // x0,becasuse read_next_frame,not read_next_x^M
        int i,j,k,molN,m;
        double center[3], tmpPos,tmpCenter,halfLbox[3],totMNC=0;
        for(k=0; k<3; k++){
                halfLbox[k] = Lbox[k]/2;
                }
        for(j=0; j<atoms; j++){
                totMNC += mass[j];
        }
        for(i=0; i<nm; i++)
        {
                molN = index[i*atoms];
                for(k=0; k<3; k++)
                {
                        tmpCenter = 0.0;  // initiate for each molecule^M
                        for(j=0; j<atoms; j++)
                        {
                                tmpPos =fr.x[molN+j][k];
                                // first unmap the atoms of big molecule into one whole molecule (use 1st atom--x0[molN][k] as absolute location)^M
                                while(tmpPos-fr.x[molN][k]> halfLbox[k])
                                        tmpPos -= Lbox[k];
                                while(tmpPos-fr.x[molN][k]<-halfLbox[k])
                                        tmpPos += Lbox[k];
                                tmpCenter += tmpPos*mass[j];
                        }
                        posc[i][k] = tmpCenter/totMNC;
                        // move molecules into box.^M
                        while(posc[i][k]>Lbox[k])
                                posc[i][k] -= Lbox[k];
                        while(posc[i][k]<0)
                                posc[i][k] += Lbox[k];
                }
        }
        return 0;
}


void getMolMultipole_z(t_trxframe fr, int atoms, int nm, int *index, double **posc, double Lbox[], double* divC, double *charge,
                                        double *dipole, double *quadrupole)
{
        int i, molN, j;
        double zmj,totC, totP, totQ;
        double A_of_electtrode = 4.176 * 4.254;
        for (i=0;i<nm;i++)
        {
                molN=index[i*atoms];
                totC=0;
                totP=0;
                totQ=0;
                for(j=0;j<atoms;j++)
                {
                        zmj=fr.x[molN+j][2]-posc[i][2];
                        zmj=pbc(zmj,Lbox[2]);
                        totC += divC[j];
                        totP += divC[j] * zmj;
                        totQ += 0.5*divC[j] * zmj *zmj;
                }
                charge[i] = totC/ A_of_electtrode;
                dipole[i] = totP/ A_of_electtrode;
                quadrupole[i] = totQ/ A_of_electtrode;
        }
}

int main(int argc, char** argv)
{
        return gmx_run_cmain(argc, argv, &main_func);
}