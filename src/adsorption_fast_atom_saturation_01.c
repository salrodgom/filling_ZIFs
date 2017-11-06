/****************************************************************************

-----------------------------------------------------------------------------
---------------     adsorption_fast_atom_saturation_01   --------------------
-----------                                                     -------------
----           (a) Univ. Pablo de Olabvide, Seville, Spain               ----
--------------------------    April 2017       -----------------------------                              
-----------------------------------------------------------------------------

  The code fills the void space of a pòrous solids with Ar atoms, following
  the crystal structure of solid Ar
    
  crash_label is an error message flag: = 0 if error ocurrs, > 0 else 
*****************************************************************************/

# include <math.h>
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>

# define MaxAtomNum   50000
# define MaxAtomNum1  20000
# define MaxAtomNum2  35000
# define MaxAtomNum3  10000
# define MaxAtomNum4  500
# define MaxAtomNum5  50
# define MaxAtomNum6  205

# define Maxline      300
# define N            10  
# define TOidealDist  1.605
# define TidealAngle  109.47122
  
FILE  *InpF;
FILE  *OutpF;
char   InpFName[Maxline], OutpFName[Maxline];
double xatom[MaxAtomNum1], yatom[MaxAtomNum1], zatom[MaxAtomNum1];
char   atomname[MaxAtomNum1][8];
char   atom_c_s[MaxAtomNum1];
double atomcharge[MaxAtomNum1];
int    CoorFType;
int    atomnumber;
double cell[N];
double Arx0[N], Ary0[N], Arz0[N];
double Arx[MaxAtomNum1], Ary[MaxAtomNum1], Arz[MaxAtomNum1];
double Arxca[MaxAtomNum1], Aryca[MaxAtomNum1], Arzca[MaxAtomNum1];
int    a_supercell, b_supercell, c_supercell;
int    minloading, maxloading;
int    removeAr[MaxAtomNum1];
double Arcell0[N], Arcell[N], Arcell_displace[4], solidcell_displace[4], tol_dist_Ar_solid;
int    Aratomnumber;
double Arxout[MaxAtomNum1], Aryout[MaxAtomNum1], Arzout[MaxAtomNum1];
double Arxtmp[MaxAtomNum1], Arytmp[MaxAtomNum1], Arztmp[MaxAtomNum1];
int    Arneighb[MaxAtomNum1][MaxAtomNum3], numbArneighb[MaxAtomNum1];
double tol_dist_Ar_Ar, dist_Ar_Ar;

int    axesi[4], icellstat[MaxAtomNum6][MaxAtomNum6][MaxAtomNum6];
int    icellframeatom[MaxAtomNum6][MaxAtomNum6][MaxAtomNum6];
double axesidelta[4], icelldframeatom[MaxAtomNum6][MaxAtomNum6][MaxAtomNum6];
double icell_1pr[MaxAtomNum6][MaxAtomNum6][MaxAtomNum6];
int    accessatom[MaxAtomNum1], numaccessatom;
int    nsorticellstat;
int  **sorticellstat;
int    tsorticellstat[MaxAtomNum][4];
int    Ar1atomnumber;
double Ar1x[MaxAtomNum3], Ar1y[MaxAtomNum3], Ar1z[MaxAtomNum3];
int  **sortArdist;

void   Initialization();
void   MainJob();
double Cryst2Cartes(int i, double cell[N], double per_atom[N]);
double Cartes2Cryst(int i, double cell[N], double per_atom[N]);
int    getatomcoor(char coordFName[Maxline], int CoorFType);
char  *grepnword(char *line1, int position_line, int nthword);
char  *grepnmword(char *line1, int position_ini, int position_final);
int    strgrepf(char line1[], char line2[],int maxline1, FILE *filename1);
int    countwords(char line1[Maxline], int count1pos, int count2pos);
double atomsdistc(double x1, double y1, double z1, double x2, double y2, double z2);
double atomsdist (double x1, double y1, double z1, double x2, double y2, double z2);
double absx(double x);

int    icellneighbcount(int ix, int iy, int iz);
int    compare_qsort( const void *pa, const void *pb);
//int    compare_qsort( const void *pa, const void *pb, const void *pc, const void *pd); 
void   icelloverlap_Ar_atom(int Ndisp_a1, int Ndisp_b1, int Ndisp_c1, int Ar_or_atom, double xatom1, double yatom1, double zatom1);

int main(argc, argv)
        int argc;
        char **argv;
{
float xx1, yy1, zz1, alfa,  xx2, yy2, zz2, xx3, yy3, zz3, xrotaxis, yrotaxis, zrotaxis;
float xrot1, yrot1, zrot1, d1, d2, d3;
int   count1, count2, count3, MainCount, crash_label; 
char  line[Maxline];
int   maxline;

     maxline = Maxline;

     crash_label = 0.0;
     
     memset(InpFName, '\0', sizeof(InpFName));
     strcpy(InpFName, argv[1]);

     if (argc>=3) 
	 {
	   memset(line, '\0', sizeof(line));
	   strcpy(line, argv[2]);
	   if ( (strcmp(line, "1") == 0) || (strcmp(line, "2") == 0) || (strcmp(line, "3") == 0) || (strcmp(line, "4") == 0) || (strcmp(line, "5") == 0) )
	   {
	     sscanf(line,"%d", &CoorFType); 	
	   }
	   else
	   {
	     CoorFType = 5;
	   }
	 }  
	 else
	 {
	   CoorFType = 5;
	 }   
     
	 //InpF = fopen(argv[1],"r");
         
     //printf(" Processing  %s CoorFType  %d\n", InpFName, CoorFType);
     
     Initialization();
     
     MainJob();
 
     printf(" thanks for using adsorption_fast_atom_saturation_01 \n");
    
}     /******************************  end of main  **************************/

/******************************************************************************
function to initiate the job; reading general data and initialising
general variables.

******************************************************************************/
void    Initialization()
{   
    int  maxline, count1, count2;
    char line[Maxline], s1[Maxline], s2[Maxline];

    maxline = Maxline;     
    
    Arx0[1] = 0.0; Ary0[1] = 0.0; Arz0[1] = 0.0; 
    Arx0[2] = 0.0; Ary0[2] = 0.5; Arz0[2] = 0.5; 
    Arx0[3] = 0.5; Ary0[3] = 0.0; Arz0[3] = 0.5;
    Arx0[4] = 0.5; Ary0[4] = 0.5; Arz0[4] = 0.0;  
    
    Arcell0[1] = 5.256; Arcell0[2] = 5.256; Arcell0[3] = 5.256;
    Arcell0[4] = 90.0;  Arcell0[5] = 90.0;  Arcell0[6] = 90.0;
    dist_Ar_Ar        = Arcell0[1]*sqrt(2.0)/2.0; 
    tol_dist_Ar_solid = 0.588 + Arcell0[1]*sqrt(2.0)/4.0;    // average atom size + Ar size
    // 47 Ar atoms are found for ZIF8_RASPA.cif using 0.725 tol_dist_Ar_solid = 0.725 + Arcell0[1]*sqrt(2.0)/4.0;
   
    //for (count1 = 0; count1 < MaxAtomNum; count1++)
    //{
    //  removeAr[count1] = 0;
    //}

    for (count1 = 0; count1 < MaxAtomNum1; count1++)
    {
      numbArneighb[count1] = -1;
      for (count2 = 0; count2 < MaxAtomNum3; count2++)
      {
        Arneighb[count1][count2] = -1;
      } 
    }

    memset(s1, '\0', sizeof(s1));
    memset(s2, '\0', sizeof(s2));
    strcpy(s1, InpFName);
    count2 = -1;
    for (count1 = strlen(s1) -1; count1 >= 0; count1--)
    {
      if (s1[count1] == '.')
      {
        count2 = count1;
        count1 = -1;
      }
    }
    if (count2 > 0)
    {
      strncpy(s2, s1, count2);
      strcpy(OutpFName, s2);
      strcat(OutpFName,"_Ar.cif");
    }
    else
    {
      strcpy(OutpFName, InpFName);
      strcat(OutpFName,"_Ar.cif");
      strcpy(s2, InpFName);
    }
    
    OutpF = fopen(OutpFName,"w");
    fprintf(OutpF, "#######################################################################\n");
    fprintf(OutpF, "# Data File\n");
    fprintf(OutpF, "# Maximun Ar loaded porous solid %s\n", InpFName);
    fprintf(OutpF, "#######################################################################\n");
    fprintf(OutpF, "data_%s_Ar\n", s2);
    fprintf(OutpF, "\n");
    fprintf(OutpF, "_audit_creation_method            'adsorption_fast_atom_saturation contact Ruiz-Salvador, Balestra or Hamad'\n");
    fprintf(OutpF, "\n");
    fprintf(OutpF, "_symmetry_space_group_name_H-M    'P1'\n");
    fprintf(OutpF, "_symmetry_Int_Tables_number       1\n");
    fprintf(OutpF, "_symmetry_cell_setting            triclinic\n");
    fprintf(OutpF, "\n");
    fprintf(OutpF, "loop_\n");
    fprintf(OutpF, "_symmetry_equiv_pos_as_xyz\n");
    fprintf(OutpF, " x,y,z\n");
    fprintf(OutpF, "\n");
    
}    /***********    end of Initialization     ************/

/******************************************************************************
function to do the Main Job
******************************************************************************/
void   MainJob()
{
    char line[Maxline];
    int    count1, count2, count3, count4, count5, count6, count7, count8, maxline;
    int    tmpAtoms[MaxAtomNum];
    double disttmp;
    double delta_a[4], delta_b[4], delta_c[4], Arcell[N];
    double value1, value2, value3;
    double shift_x, shift_y, shift_z, delta_shift;
    double Ar_box_disp[4], Arboxcentre[4], solidcentre[4];
    double cryst_atom[4], tmpcell[N];
    int    disp_count1, disp_count2, disp_count3;
    int    Ndisp_a, Ndisp_b, Ndisp_c;
    double small_shift_a, small_shift_b, small_shift_c;
    double x_tmp, y_tmp, z_tmp;
    double disttmp1;
    
    printf(" ms 1 Processing  %s CoorFType  %d Output data to %s\n", InpFName, CoorFType, OutpFName);
   
    count1 = getatomcoor(InpFName, CoorFType);
    
    fprintf(OutpF, "_cell_length_a                    %.5f\n", cell[1]);
    fprintf(OutpF, "_cell_length_b                    %.5f\n", cell[2]);
    fprintf(OutpF, "_cell_length_c                    %.5f\n", cell[3]);
    fprintf(OutpF, "_cell_angle_alpha                 %.5f\n", cell[4]);
    fprintf(OutpF, "_cell_angle_beta                  %.5f\n", cell[5]);
    fprintf(OutpF, "_cell_angle_gamma                 %.5f\n", cell[6]);
    fprintf(OutpF, "\n");
    fprintf(OutpF, "loop_\n");
    fprintf(OutpF, "_atom_site_label\n");
    fprintf(OutpF, "_atom_site_type_symbol\n");
    fprintf(OutpF, "_atom_site_fract_x\n");
    fprintf(OutpF, "_atom_site_fract_y\n");
    fprintf(OutpF, "_atom_site_fract_z\n");
    fprintf(OutpF, "_atom_site_charge\n");

    for (count1 = 0; count1 < atomnumber; count1++)
    {
      fprintf(OutpF, " %s%d %s %.5f %.5f %.5f 0.00\n", atomname[count1], count1+1, atomname[count1], xatom[count1], yatom[count1], zatom[count1]);
    }
    
    for (count1 = 0; count1 < atomnumber; count1++)
    {
      if (xatom[count1] < 0.0)
      {
        xatom[count1] = xatom[count1] + 1.0;
      }
      else
      {
      	if (xatom[count1] >= 1.0)
        {
          xatom[count1] = xatom[count1] - 1.0;
        }
      }
      if (yatom[count1] < 0.0)
      {
        yatom[count1] = yatom[count1] + 1.0;
      }
      else
      {
      	if (yatom[count1] >= 1.0)
        {
          yatom[count1] = yatom[count1] - 1.0;
        }
      }
      if (zatom[count1] < 0.0)
      {
        zatom[count1] = zatom[count1] + 1.0;
      }
      else
      {
      	if (zatom[count1] >= 1.0)
        {
          zatom[count1] = zatom[count1] - 1.0;
        }
      }
    }

    //for (count1 = 0; count1 <= atomnumber; count1++)
    //{
    //  accessatom[count1] = -1;
    //}
    //numaccessatom = 0;
    
    for (count1 = 1; count1 <= 3; count1++)
    {
      if (cell[count1] > 100)
      {
        axesi[count1]      = 200;
        axesidelta[count1] = cell[count1]/200;
      }
      else
      {
        axesi[count1]      = round(cell[count1]/0.5);
        axesidelta[count1] = cell[count1]/axesi[count1];
      }
    }
    //printf(" ms 1b %d %f axesi-x %d %f axesi-y %d %f axesi-z\n", axesi[1], axesidelta[1], axesi[2], axesidelta[2], axesi[3], axesidelta[3]); 


// OJO OJO aaaa aqui empieza lo ya bueno, cuando aun no se habian encontrado las posiciones de los Ar1 

    // computing delta_a
    if (cell[6] > 90.0)
    {
      delta_a[2] = -cell[2]*cos(cell[6]*M_PI/180.0);
      delta_b[2] = -cell[1]*cos(cell[6]*M_PI/180.0);
    }
    else
    {
      delta_a[2] = cell[2]*cos(cell[6]*M_PI/180.0);
      delta_b[2] = cell[1]*cos(cell[6]*M_PI/180.0);
    }
    if (cell[5] > 90.0)
    {
      delta_a[3] = -cell[3]*cos(cell[5]*M_PI/180.0);
      delta_c[3] = -cell[1]*cos(cell[5]*M_PI/180.0);
    }
    else
    {
      delta_a[3] = cell[3]*cos(cell[5]*M_PI/180.0);
      delta_c[3] = cell[1]*cos(cell[5]*M_PI/180.0);
    }
    if (cell[4] > 90.0)
    {
      delta_b[2] = -cell[3]*cos(cell[5]*M_PI/180.0);
      delta_c[2] = -cell[2]*cos(cell[5]*M_PI/180.0);
    }
    else
    {
      delta_b[2] = cell[3]*cos(cell[5]*M_PI/180.0);
      delta_c[2] = cell[2]*cos(cell[5]*M_PI/180.0);
    }
    delta_a[1] = delta_a[2] + delta_a[3];
    delta_b[1] = delta_b[2] + delta_b[3];
    delta_c[1] = delta_c[2] + delta_c[3];
  
    //printf("\n");
    //printf(" ms2 cell parameters %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);
    //printf(" ms2a delta_a %.5f %.5f %.5f\n", delta_a[1], delta_b[1], delta_c[1]);
    //printf(" ms2b delta_a %.5f %.5f %.5f\n", delta_a[2], delta_b[2], delta_c[2]);
    //printf(" ms2c delta_a %.5f %.5f %.5f\n", delta_a[3], delta_b[3], delta_c[3]);

    
    value1 = cell[1] + delta_a[1];
    a_supercell = (value1 / Arcell0[1]) + 2;
    value1 = cell[2] + delta_b[1];
    b_supercell = (value1 / Arcell0[2]) + 2;
    value1 = cell[3] + delta_c[1];
    c_supercell = (value1 / Arcell0[3]) + 2;
    Arcell[1] = Arcell0[1] * a_supercell; Arcell[2] = Arcell0[2] * b_supercell; Arcell[3] = Arcell0[3] * c_supercell;
    Arcell[4] = Arcell0[4]; Arcell[5] = Arcell0[5]; Arcell[6] = Arcell0[6]; 
  
    // print Arcell and delta_cell related information for checking
    //printf(" ms3a Ar0_cell = %.5f supercell: %.5f %.5f %.5f %.5f %.5f %.5f\n", Arcell0[1], Arcell[1], Arcell[2], Arcell[3], Arcell[4], Arcell[5], Arcell[6]);
    //printf(" ms3b Ar supercell factors: %d %d %d\n", a_supercell, b_supercell, c_supercell); 
    //printf("\n");
 
    //finding displacement vector of the Ar box centre to the solid centre
    //Ar_box_disp[4], Arboxcentre[4], solidcentre[4]; 
    cryst_atom[1]  = 0.5; 
    cryst_atom[2]  = 0.5;
    cryst_atom[3]  = 0.5;
    solidcentre[1] = Cryst2Cartes(1, cell, cryst_atom);
    solidcentre[2] = Cryst2Cartes(2, cell, cryst_atom);
    solidcentre[3] = Cryst2Cartes(3, cell, cryst_atom);
    tmpcell[1] = cell[1]; tmpcell[2] = cell[2]; tmpcell[3] = cell[3]; 
    tmpcell[4] = cell[4]; tmpcell[5] = cell[5]; tmpcell[6] = cell[6];
    cell[1] = Arcell[1];  cell[2] = Arcell[2];  cell[3] = Arcell[3];
    cell[4] = Arcell[4];  cell[5] = Arcell[5];  cell[6] = Arcell[6];
    cryst_atom[1]  = 0.5; 
    cryst_atom[2]  = 0.5;
    cryst_atom[3]  = 0.5;
    Arboxcentre[1] = Cryst2Cartes(1, cell, cryst_atom) - Arcell0[1]/2.0;
    Arboxcentre[2] = Cryst2Cartes(2, cell, cryst_atom) - Arcell0[2]/2.0;
    Arboxcentre[3] = Cryst2Cartes(3, cell, cryst_atom) - Arcell0[3]/2.0;
    Ar_box_disp[1] = solidcentre[1] - Arboxcentre[1];
    Ar_box_disp[2] = solidcentre[2] - Arboxcentre[2];
    Ar_box_disp[3] = solidcentre[3] - Arboxcentre[3];

    //printf(" ms4 Ar_box_centre %.5f %.5f %.5f\n", Arboxcentre[1], Arboxcentre[2], Arboxcentre[3]);
    //printf(" ms4 solidcentre   %.5f %.5f %.5f\n", solidcentre[1], solidcentre[2], solidcentre[3]);
    //printf(" ms4 Ar_box_disp   %.5f %.5f %.5f\n", Ar_box_disp[1], Ar_box_disp[2], Ar_box_disp[3]);
    //printf(" ms4a Ar0_1 %.5f %.5f %.5f\n", Arx0[1], Ary0[1], Arz0[1]);
    //printf(" ms4a Ar0_2 %.5f %.5f %.5f\n", Arx0[2], Ary0[2], Arz0[2]);
    //printf(" ms4a Ar0_3 %.5f %.5f %.5f\n", Arx0[3], Ary0[3], Arz0[3]);
    //printf(" ms4a Ar0_4 %.5f %.5f %.5f\n", Arx0[4], Ary0[4], Arz0[4]);

    //printf("\n");

    // OJO AAAA aqui empieza lo ya bueno        

    //creating Ar box, displace Ar atoms 
    Aratomnumber = 0;
    for (count1 = 1; count1 <= a_supercell; count1++)
    {
      for (count2 = 1; count2 <= b_supercell; count2++)
      {
        for (count3 = 1; count3 <= c_supercell; count3++)
        {
          for (count4 = 1; count4 <= 4; count4++)
          {	
            Aratomnumber++;
            //cryst_atom[1]       = Arx0[count4]/a_supercell + count1 - 1; 
            //cryst_atom[2]       = Ary0[count4]/b_supercell + count2 - 1;
            //cryst_atom[3]       = Arz0[count4]/c_supercell + count3 - 1;
            cryst_atom[1]       = (Arx0[count4] + count1 - 1)/a_supercell; 
            cryst_atom[2]       = (Ary0[count4] + count2 - 1)/b_supercell;
            cryst_atom[3]       = (Arz0[count4] + count3 - 1)/c_supercell;
            Arxca[Aratomnumber] = Cryst2Cartes(1, cell, cryst_atom) + Ar_box_disp[1];
            Aryca[Aratomnumber] = Cryst2Cartes(2, cell, cryst_atom) + Ar_box_disp[2];
            Arzca[Aratomnumber] = Cryst2Cartes(3, cell, cryst_atom) + Ar_box_disp[3];
          } 
        }      
      }
    }
    // OJO only for checking debuging purpose, comment on standard running
    //for (count1 = 1; count1 <= Aratomnumber; count1++)
    //{
    //  printf(" ms5 %d\t Ar cartes %.5f %.5f %.5f\n", count1, Arxca[count1], Aryca[count1], Arzca[count1]);
    //  printf(" ms5 %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);
    //}
    //printf("\n");

    cell[1] = tmpcell[1]; cell[2] = tmpcell[2]; cell[3] = tmpcell[3];
    cell[4] = tmpcell[4]; cell[5] = tmpcell[5]; cell[6] = tmpcell[6];
    // first scan to delete part of Ar atoms, those out and overlaping the solid and
    // recovering back those at distances Arcell0/2 below of already labeled good Ar atoms
    for (count1 = 1; count1 <= Aratomnumber; count1++)
    {
      removeAr[count1] = 0;
      Arxout[count1] = Arx[count1]; Aryout[count1] = Ary[count1]; Arzout[count1] = Arz[count1];
    }
    for (count1 = 1; count1 <= Aratomnumber; count1++)
    {
      cryst_atom[1] = Arxca[count1];
      cryst_atom[2] = Aryca[count1];
      cryst_atom[3] = Arzca[count1];
      Arx[count1]   = Cartes2Cryst(1, cell, cryst_atom);
      Ary[count1]   = Cartes2Cryst(2, cell, cryst_atom);
      Arz[count1]   = Cartes2Cryst(3, cell, cryst_atom);
      if ( ((Arx[count1]>=0.0) && (Arx[count1]<=1.0)) && ((Ary[count1]>=0.0) && (Ary[count1]<=1.0)) && ((Arz[count1]>=0.0) && (Arz[count1]<=1.0)) )
      {
        //Ar atom in the solid and now check overlap
        for (count2 = 0; count2 < atomnumber; count2++)
        {
          disttmp = atomsdist(xatom[count2], yatom[count2], zatom[count2], Arx[count1], Ary[count1], Arz[count1]);
          if (disttmp < tol_dist_Ar_solid)
          {
            removeAr[count1] = 1;
            count2 = atomnumber + 1;
          }
        }
      }
      else
      {
          //fuera del solido
          removeAr[count1] = 2;
      }
    }
    count5 = 0;
    for (count1 = 1; count1 <= Aratomnumber; count1++)
    {
      if (removeAr[count1] == 0)
      {
        count5++;
        Arxtmp[count5] = Arxca[count1]; Arytmp[count5] = Aryca[count1]; Arztmp[count5] = Arzca[count1];
      }
    }
    count4 = count5;  
    //printf(" ms5a %d\t Ar after first cut (out of the solid and overlap solid atoms, cutoff %.5f)\n", count5, tol_dist_Ar_solid);
    //printf(" ms5a %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);
    //for (count1 = 1; count1 <= count4; count1++)
    //{
    //  printf(" ms5a %d\t Ar cartes %.5f %.5f %.5f\n", count1, Arxtmp[count1], Arytmp[count1], Arztmp[count1]);
    //}
    //printf("\n");
     
    tol_dist_Ar_Ar = dist_Ar_Ar*1.05;
    for (count1 = 1; count1 <= Aratomnumber; count1++)
    {
      if (removeAr[count1] == 1)
      {
        for (count2 = 1; count2 <= count4; count2++) 
        {
          disttmp = atomsdistc(Arxtmp[count2], Arytmp[count2], Arztmp[count2], Arxca[count1], Aryca[count1], Arzca[count1]);
          if (disttmp <= tol_dist_Ar_Ar) 
          // if (disttmp <= Arcell0[1]*1.05) aaa try this if does not work
          {
            removeAr[count1] = 0;
            count2 = count4 + 1;
            count5++;
            Arxtmp[count5] = Arxca[count1]; Arytmp[count5] = Aryca[count1]; Arztmp[count5] = Arzca[count1];
          }
        }
      }
    }
    //count5 = 0; 
    //for (count1 = 1; count1 <= Aratomnumber; count1++)
    //{
    //  if (removeAr[count1] == 0)
    //  {
    //    count5++;
    //    Arxtmp[count5] = Arxca[count1]; Arytmp[count5] = Aryca[count1]; Arztmp[count5] = Arzca[count1];
    //  }
    //}
    Aratomnumber = count5;
    //printf(" ms5b %d\t Ar after addition after first cut and without deleting duplicated atoms\n", count5);
    //printf(" ms5b %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);
    for (count1 = 1; count1 <= count5; count1++)
    { 
      Arxca[count1] = Arxtmp[count1]; Aryca[count1] = Arytmp[count1]; Arzca[count1] = Arztmp[count1];
      //printf(" ms5b %d\t Ar cartes %.5f %.5f %.5f\n", count1, Arxca[count1], Aryca[count1], Arzca[count1]);
    }
    //printf("\n");

    // cleaning duplicate Ar atoms
    //  NOT NEEDED NOW
    //
    //for (count1 = 1; count1 <= Aratomnumber; count1++)
    //{
    //  removeAr[count1] = 0;
    //} 
    //for (count1 = 2; count1 <= Aratomnumber; count1++)
    //{
    //  for (count2 = 1; count2 <count1; count2++)
    //  {
    //    disttmp = atomsdistc(Arxtmp[count2], Arytmp[count2], Arztmp[count2], Arxtmp[count1], Arytmp[count1], Arztmp[count1]);
    //    //if (disttmp < Arcell0[1]*0.75*sqrt(2.0)/2.0)
    //    if (disttmp < 0.5)  
    //    {
    //      removeAr[count1] = 1;
    //      printf(" ms5z Aratom %d %.5f %.5f %.5f ", count1, Arxtmp[count1], Arytmp[count1], Arztmp[count1]);
    //      printf("Artmp %d %.5f %.5f %.5f\n ", count2, Arxtmp[count2], Arytmp[count2], Arztmp[count2]);
    //      count2 = count1 + 1;
    //    }
    //  } 
    //}
    //count5 = 0;
    //for (count1 = 1; count1 < Aratomnumber; count1++)
    //{
    //  if (removeAr[count1] == 0)
    //  {
    //    count5++;
    //    Arxca[count5] = Arxtmp[count1]; Aryca[count5] = Arytmp[count1]; Arzca[count5] = Arztmp[count1];
    //  }
    //}
    //Aratomnumber = count5;
    //printf(" ms5d %d\t Ar after addition after first cut and WITH    deleting duplicated atoms\n", count5);
    //printf(" ms5d %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);

    // qsort of Ar atoms according to their distance with the center of the solid (and displaced Ar template)
    sortArdist = malloc(Aratomnumber * sizeof(int*));
    for (count1 = 0; count1 < Aratomnumber; count1++)
    {
      count2                = round(10000.0*atomsdistc(solidcentre[1], solidcentre[2], solidcentre[3], Arxca[count1 + 1], Aryca[count1 + 1], Arzca[count1 + 1]));
      sortArdist[count1]    = malloc(2 * sizeof(int));
      sortArdist[count1][0] = count2;
      sortArdist[count1][1] = count1 + 1;
    } 
                
    qsort(sortArdist, Aratomnumber, sizeof sortArdist[0], compare_qsort);
    
    for (count1 = 0; count1 < Aratomnumber; count1++)
    {
      count2 = sortArdist[count1][1];
      Arxtmp[count1+1] = Arxca[count2]; Arytmp[count1+1] = Aryca[count2]; Arztmp[count1+1] = Arzca[count2];
    }
    for (count1 = 1; count1 <= Aratomnumber; count1++)
    {
      Arxca[count1] = Arxtmp[count1]; Aryca[count1] = Arytmp[count1]; Arzca[count1] = Arztmp[count1];
      //printf(" ms6 %d\t Ar cartes %.5f %.5f %.5f\n", count1, Arxca[count1], Aryca[count1], Arzca[count1]);
    }

    // finding the neighbours of the solid atoms for each Ar atom, tolerance distance = 1.0*Arcell0
    //tol_dist_Ar_Ar = Arcell0[1]*1.15*sqrt(3.0);
    //tol_dist_Ar_Ar = dist_Ar_Ar;
    tol_dist_Ar_Ar = Arcell0[1];
    for (count1 = 1; count1 <= Aratomnumber; count1++)
    {
      cryst_atom[1] = Arxca[count1];
      cryst_atom[2] = Aryca[count1];
      cryst_atom[3] = Arzca[count1];
      value1        = Cartes2Cryst(1, cell, cryst_atom);
      value2        = Cartes2Cryst(2, cell, cryst_atom);
      value3        = Cartes2Cryst(3, cell, cryst_atom);
      count3 = 0; // counts numbArneighb
      for (count2 = 0; count2 < atomnumber; count2++)
      {
        disttmp = atomsdist(xatom[count2], yatom[count2], zatom[count2], value1, value2, value3);
        if (disttmp <= tol_dist_Ar_Ar) 
        {
          count3++;
          Arneighb[count1][count3] = count2;
        }
      }
      numbArneighb[count1] = count3;
      //printf(" ms7 Ar %d\t numbArneighb %d\n", count1, numbArneighb[count1]); 
    }

    // verify in cell Ar atoms and eliminate overlaping atoms with the solid
    // displace the Ar solid on aroun 1/2 of its cell parameter
    Ndisp_a =  4; Ndisp_b =  4; Ndisp_c =  4;
    //Ndisp_a = 10; Ndisp_b = 10; Ndisp_c = 10;
    //Ndisp_a = 1; Ndisp_b = 1; Ndisp_c = 1; //OJO comentar es solo para test
    small_shift_a = Arcell0[1]*sqrt(2.0)/2.0 / (Ndisp_a*1.0); 
    small_shift_b = Arcell0[2]*sqrt(2.0)/2.0 / (Ndisp_b*1.0); 
    small_shift_c = Arcell0[3]*sqrt(2.0)/2.0 / (Ndisp_c*1.0);
    //small_shift_a = Arcell0[1] / Ndisp_a; small_shift_b = Arcell0[2] / Ndisp_b; small_shift_c = Arcell0[3] / Ndisp_c;
    minloading = 10000; maxloading = -1; count5 = -1;
    cell[1] = tmpcell[1]; cell[2] = tmpcell[2]; cell[3] = tmpcell[3];
    cell[4] = tmpcell[4]; cell[5] = tmpcell[5]; cell[6] = tmpcell[6];
    tol_dist_Ar_Ar = dist_Ar_Ar*0.95;
    for (disp_count1 = -Ndisp_a; disp_count1 <= Ndisp_a; disp_count1++)
    {
      //printf(" ms8 Ndisp_a %d minloading %d maxloading %d Aratomnumber %d loading %d tol_dist_Ar_solid %.5f\n", disp_count1, minloading, maxloading, Aratomnumber, count5, tol_dist_Ar_solid); 
      for (disp_count2 = -Ndisp_b; disp_count2 <= Ndisp_b; disp_count2++)
      {
        //printf(" ms9 Ndisp_b %d Aratomnumber %d loading %d\n", disp_count2, Aratomnumber, count5); 
        for (disp_count3 = -Ndisp_c; disp_count3 <= Ndisp_c; disp_count3++)
        {
          //printf(" ms10 Ndisp_c %d Aratomnumber %d loading %d\n", disp_count3, Aratomnumber, count5);
          //printf(" ms10 %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);
          count4 = 0; count6 = 0; count7 = 0; count8 = 0;
          for (count1 = 1; count1 <= Aratomnumber; count1++)
          {
            removeAr[count1] = 0;
            cryst_atom[1] = Arxca[count1] + 0.5*small_shift_a*disp_count1;
            cryst_atom[2] = Aryca[count1] + 0.5*small_shift_b*disp_count2;
            cryst_atom[3] = Arzca[count1] + 0.5*small_shift_c*disp_count3;
            Arx[count1]   = Cartes2Cryst(1, cell, cryst_atom);
            Ary[count1]   = Cartes2Cryst(2, cell, cryst_atom);
            Arz[count1]   = Cartes2Cryst(3, cell, cryst_atom);
            //printf(" ms11a Ar %3d %3d %3d solid_displacement %.5f %.5f %.5f\n", disp_count1, disp_count2, disp_count3, Arx[count1], Ary[count1], Arz[count1]);       
          }
          for (count1 = 1; count1 <= Aratomnumber; count1++)
          {
            if (removeAr[count1] == 0)
            { 
              //printf(" ms11a %d\t Ar cryst ?? solid %.5f %.5f %.5f\n", count1, Arx[count1], Ary[count1], Arz[count1]);
              if ( ((Arx[count1]>=0.0) && (Arx[count1]<=1.0)) && ((Ary[count1]>=0.0) && (Ary[count1]<=1.0)) && ((Arz[count1]>=0.0) && (Arz[count1]<=1.0)) )
              {
                count4++;
                //printf(" ms11b %d\t Ar cryst IN solid %.5f %.5f %.5f\n", count1, Arx[count1], Ary[count1], Arz[count1]);
                //Ar atom in the solid and now check overlap
                for (count2 = 1; count2 <= numbArneighb[count1]; count2++)
                {
                  count3 = Arneighb[count1][count2];
                  disttmp = atomsdist(xatom[count3], yatom[count3], zatom[count3], Arx[count1], Ary[count1], Arz[count1]);     	
                  if (disttmp < tol_dist_Ar_solid)
                  {
                    removeAr[count1] = 1;
                    count2 = atomnumber + 1;
                    count6++;
                  }
                }
              }
              else
              {
                //fuera del solido
                removeAr[count1] = 1;
                count8++;
              }
            }
            if (removeAr[count1] == 0)
            {
              //hacer un ciclo primero para eliminar y marcar como removeAr[count1] = 1 todos los Ar unidos al corriente
              for (count2 = count1 + 1; count2 <= Aratomnumber; count2++)
              {
                if (removeAr[count2] == 0)
                {
                  disttmp = atomsdist(Arx[count2], Ary[count2], Arz[count2], Arx[count1], Ary[count1], Arz[count1]);
                  if (disttmp < tol_dist_Ar_Ar) 
                  //if (disttmp < 0.5)
                  {
                    removeAr[count2] = 1;
                    count7++;
                    //printf(" ms11b Ar %3d number %3d removed %3d %3d %3d solid_displacement %.5f %.5f %.5f coord as dist Ar-Ar %.5f\n", count2, count7, disp_count1, disp_count2, disp_count3, Arxtmp[count5], Arytmp[count5], Arztmp[count5], disttmp);
                  }
                }
              }
            }
          }
          count5 = 0;
          for (count1 = 1; count1 < Aratomnumber; count1++)
          {
            if (removeAr[count1] == 0)
            {
              count5++;
              Arxtmp[count5] = Arx[count1]; Arytmp[count5] = Ary[count1]; Arztmp[count5] = Arz[count1];
              //printf(" ms11c Ar %3d %3d %3d solid_displacement %.5f %.5f %.5f\n", disp_count1, disp_count2, disp_count3, Arxtmp[count5], Arytmp[count5], Arztmp[count5]);
            }
          }
          if (count5 <= minloading)
          {
            minloading = count5;
          }
          if (count5 > maxloading)
          {
            maxloading = count5;
            for (count1 = 1; count1 <= maxloading; count1++)
            {
              Arxout[count1] = Arxtmp[count1]; Aryout[count1] = Arytmp[count1]; Arzout[count1] = Arztmp[count1];
              //printf(" ms11d Ar %3d %3d %3d solid_displacement %.5f %.5f %.5f maxloading %3d\n", disp_count1, disp_count2, disp_count3, Arxout[count1], Aryout[count1], Arzout[count1], maxloading);
            }
          }
          //printf(" ms12a %d out_solid_cell %d overlaped_solid_atoms %d duplicated %d in_solid_cell %d current_loading %d Ar_number\n", count8, count6, count7, count4, count5, Aratomnumber);
        }
      }
    }

    // only for checking debuging purpose, comment on standard running
    //for (count1 = 1; count1 <= maxloading; count1++)
    //{
    //  printf(" ms13 Ar %d\t Loading Min-Max %d-%d Coord cryst %.5f %.5f %.5f\n", count1, minloading, maxloading, Arxout[count1], Aryout[count1], Arzout[count1]);
    //}
    
    count3 = atomnumber;
    for (count1 = 1; count1 <= maxloading; count1++)
    {
      count3++;
      fprintf(OutpF, " Ar%d Ar %.5f %.5f %.5f 0.00\n", count3, Arxout[count1], Aryout[count1], Arzout[count1]);
    }
    for (count1 = 1; count1 <= maxloading; count1++)
    {
      Arx[count1] = Arxout[count1]; Ary[count1] = Aryout[count1]; Arz[count1] = Arzout[count1];
    }
    Aratomnumber = maxloading;
  // OJO AAAA aqui termina lo ya bueno    

//xxxaaa
    // creating and first update of icellstat
    // icellstat = -3 not scanned yet, -2 overlap framework, -1 overlap Ar, else number of neighbour free icells
    // max value for icellstat = 26
    for (count1 = 1; count1 <= axesi[1]; count1++)
    {
      for (count2 = 1; count2 <= axesi[2]; count2++)
      {
        for (count3 = 1; count3 <= axesi[3]; count3++)
        {
          icellstat[count1][count2][count3] = -3;
          icellframeatom[count1][count2][count3] = -1;
          icelldframeatom[count1][count2][count3] = 10000.0;
          icell_1pr[count1][count2][count3] = 0.0;
        }
      }
    }

    // looking for icell overlap with the framework atoms
    Ndisp_a = (tol_dist_Ar_solid*1.05) / axesidelta[1];
    Ndisp_b = (tol_dist_Ar_solid*1.05) / axesidelta[2];
    Ndisp_c = (tol_dist_Ar_solid*1.05) / axesidelta[3];
    //printf(" ms14 icell_range for icell-atom scan %3d %3d %3d\n", Ndisp_a, Ndisp_b, Ndisp_c);
    for (count4 = 0; count4 < atomnumber; count4++)
    {
      icelloverlap_Ar_atom(Ndisp_a, Ndisp_b, Ndisp_c, 1, xatom[count4], yatom[count4], zatom[count4]);
    }

    // OJO only for test, delete on standard running
    for (count1 = 1; count1 <= axesi[1]; count1++)
    {
      x_tmp = count1*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
      for (count2 = 1; count2 <= axesi[2]; count2++)
      {
        y_tmp = count2*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
        for (count3 = 1; count3 <= axesi[3]; count3++)
        {
          z_tmp = count3*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
          if (icellstat[count1][count2][count3] == -3)
          {
            //printf(" ms15 %d cellstat after cleaning framework-environment of icell %2d %2d %2d coord %.5f %.5f %.5f\n", icellstat[count1][count2][count3], count1, count2, count3, x_tmp, y_tmp, z_tmp);
          }
        }
      }
    } 

    // looking for icell overlap with the already incorporated Ar atoms
    Ndisp_a = (tol_dist_Ar_Ar*1.05) / axesidelta[1];
    Ndisp_b = (tol_dist_Ar_Ar*1.05) / axesidelta[2];
    Ndisp_c = (tol_dist_Ar_Ar*1.05) / axesidelta[3];
    //printf(" ms16 icell_range for icell-Ar scan %3d %3d %3d\n", Ndisp_a, Ndisp_b, Ndisp_c);
    for (count4 = 1; count4 <= Aratomnumber; count4++)
    {
      icelloverlap_Ar_atom(Ndisp_a, Ndisp_b, Ndisp_c, 2, Arx[count4], Ary[count4], Arz[count4]);
    }

    // looking for the new Ar atoms to add, now qsort according the shortest distance to framework 
    // printf ONLY the rest good !!!i  -  OJO only for test, delete on standard running
    nsorticellstat = -1;
    for (count1 = 1; count1 <= axesi[1]; count1++)
    {
      x_tmp = count1*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
      for (count2 = 1; count2 <= axesi[2]; count2++)
      {
        y_tmp = count2*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
        for (count3 = 1; count3 <= axesi[3]; count3++)
        {
          z_tmp = count3*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
          if (icellstat[count1][count2][count3] == -3)
          {
            disttmp1 = 10000.0;
            for (count4 = 0; count4 < atomnumber; count4++)
            {
              disttmp = atomsdist(xatom[count4], yatom[count4], zatom[count4], x_tmp, y_tmp, z_tmp);
              if (disttmp < disttmp1)
              {
                disttmp1 = disttmp;
                //icellstat[count1][count2][count3] = -2;
                //icelldframeatom[count1][count2][count3] = disttmp;
                //icellframeatom[count1][count2][count3]  = count4;
              }   
            }
            nsorticellstat++;
            tsorticellstat[nsorticellstat][0] = round(disttmp1*10000.0);
            tsorticellstat[nsorticellstat][1] = count1;
            tsorticellstat[nsorticellstat][2] = count2;
            tsorticellstat[nsorticellstat][3] = count3;
            //tsorticellstat[nsorticellstat][0] = icellstat[count1][count2][count3];
            //tsorticellstat[nsorticellstat][1] = round(icelldframeatom[count1][count2][count3]*10000.0);
            //tsorticellstat[nsorticellstat][2] = count1;
            //tsorticellstat[nsorticellstat][3] = count2;
            //tsorticellstat[nsorticellstat][4] = count3;
            //printf(" ms17 %d cellstat after cleaning Ar-environment of icell %2d %2d %2d coord %.5f %.5f %.5f nsorticellstat %4d dist_framework %.5f\n", icellstat[count1][count2][count3], count1, count2, count3, x_tmp, y_tmp, z_tmp, nsorticellstat, disttmp1);
          }
          else
          {
            //printf(" ms18 %d cellstat after cleaning Ar-environment of icell %2d %2d %2d coord %.5f %.5f %.5f\n", icellstat[count1][count2][count3], count1, count2, count3, x_tmp, y_tmp, z_tmp);
          }
        }
      }
    }
    nsorticellstat++;
    sorticellstat = malloc(nsorticellstat * sizeof(int*));
    for (count1 = 0; count1 < nsorticellstat; count1++)
    {
      sorticellstat[count1]    = malloc(4 * sizeof(int));
      sorticellstat[count1][0] = tsorticellstat[count1][0];
      sorticellstat[count1][1] = tsorticellstat[count1][1];
      sorticellstat[count1][2] = tsorticellstat[count1][2];
      sorticellstat[count1][3] = tsorticellstat[count1][3];
    }
    
    qsort(sorticellstat, nsorticellstat, sizeof sorticellstat[0], compare_qsort);

    for (count1 = 0; count1 < nsorticellstat; count1++)
    {
      if (sorticellstat[count1][0] < tol_dist_Ar_solid*10000.0)
      {
        removeAr[count1] = 1;
      }
      else
      {
        removeAr[count1] = 0;
      }
      Arx[count1]      = sorticellstat[count1][1]*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
      Ary[count1]      = sorticellstat[count1][2]*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
      Arz[count1]      = sorticellstat[count1][3]*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
      //printf(" ms19 qsort_icellstat %6d idist(x10000) icell %3d %3d %3d coord %.5f %.5f %.5f removeAr %d\n", sorticellstat[count1][0], sorticellstat[count1][1], sorticellstat[count1][2], sorticellstat[count1][3], Arx[count1], Ary[count1], Arz[count1], removeAr[count1]);
    }

    // looking for the new Ar atoms to add 
    count4 = 0;
    for (count1 = 0; count1 < nsorticellstat; count1++)
    {
      if (removeAr[count1] == 0)
      {
        for (count2 = count1 + 1; count2 < nsorticellstat; count2++)
        {
          if (removeAr[count2] == 0)
          {
            disttmp = atomsdist(Arx[count2], Ary[count2], Arz[count2], Arx[count1], Ary[count1], Arz[count1]);
            if (disttmp < tol_dist_Ar_Ar)
            {
              removeAr[count2] = 1;
              //printf(" ms11b Ar %3d number %3d removed %3d %3d %3d solid_displacement %.5f %.5f %.5f coord as dist Ar-Ar %.5f\n", count2, count7, disp_count1, disp_count2, disp_count3, Arxtmp[count5], Arytmp[count5], Arztmp[count5], disttmp);
            }
          }
        }
      }
      if (removeAr[count1] == 0)
      {
        count4++;
      }
    } 
 
    // only for checking debuging purpose, comment on standard running
    //for (count1 = 0; count1 < nsorticellstat; count1++)
    //{
    //  if (removeAr[count1] == 0)
    //  {
    //    printf(" ms20 Na Ar %4d Coord cryst %.5f %.5f %.5f\n", count1, Arx[count1], Ary[count1], Arz[count1]);
    //  }
    //}

    count3 = atomnumber + Aratomnumber;
    for (count1 = 0; count1 < nsorticellstat; count1++)
    {
      if (removeAr[count1] == 0)
      {
        count3++;
        fprintf(OutpF, " Ar%d Ar %.5f %.5f %.5f 0.00\n", count3, Arx[count1], Ary[count1], Arz[count1]);
      }
    }
   
    printf("\n");
    printf("Ar atoms %d\n", count3 - atomnumber); 


/*
    for (count1 = 1; count1 <= axesi[1]; count1++)
    {
      x_tmp = count1*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
      //printf(" ms 1cx %d count1 %f x_tmp\n", count1, x_tmp);
      for (count2 = 1; count2 <= axesi[2]; count2++)
      {
        y_tmp = count2*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
        //printf(" ms 1cy %d count2 %f y_tmp\n", count2, y_tmp);
        for (count3 = 1; count3 <= axesi[3]; count3++)
        {
          z_tmp = count3*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
          if (count3 == 1)
          {
            if ( (count1 == 1) && (count2 == 1) )
            {
              disttmp = 100000.0;
            }
            else
            {
              if (count1 == 1)
              {
                count5  = icellframeatom[count1][count2-1][count3];
                disttmp = atomsdist(xatom[count5], yatom[count5], zatom[count5], x_tmp, y_tmp, z_tmp);
              }
              else
              {
                count5  = icellframeatom[count1-1][count2][count3];
                disttmp = atomsdist(xatom[count5], yatom[count5], zatom[count5], x_tmp, y_tmp, z_tmp);
              }
            }
          }
          else
          {
            count5  = icellframeatom[count1][count2][count3-1];
            disttmp = atomsdist(xatom[count5], yatom[count5], zatom[count5], x_tmp, y_tmp, z_tmp);
          }
          //printf(" ms 1cz scanning %4d %4d %4d icell %f %f %f\n", count1, count2, count3, x_tmp, y_tmp, z_tmp);
          if (disttmp < tol_dist_Ar_solid)
          {
            icellstat[count1][count2][count3] = -2;
            icelldframeatom[count1][count2][count3] = disttmp;
            icellframeatom[count1][count2][count3]  = count5; 
            //if (disttmp >= tol_dist_Ar_solid-0.5)
            //{
            //  if (count5 != accessatom[numaccessatom])
            //  {
            //    for (count6 = 1; count6 <= numaccessatom; count6++)
            //    {
            //      if (count5 == accessatom[count6])
            //      {
            //        count6 = numaccessatom + 10;
            //      }
            //    }
            //    if (count6<numaccessatom + 5)
            //    {
            //      numaccessatom++;
            //      accessatom[numaccessatom] = count5;
            //    }
            //  }
            //}      
          }
          else
          {
            for (count4 = 0; count4 < atomnumber; count4++)
            {  
              disttmp = atomsdist(xatom[count4], yatom[count4], zatom[count4], x_tmp, y_tmp, z_tmp);
              if (disttmp < tol_dist_Ar_solid)
              {
                icellstat[count1][count2][count3] = -2;
                icelldframeatom[count1][count2][count3] = disttmp;
                icellframeatom[count1][count2][count3]  = count4;
                //xxx accessatom se usa ???
                //if (count4 != accessatom[numaccessatom])
                //{
                //  if (disttmp >= tol_dist_Ar_solid-0.5)
                //  {
                //    for (count6 = 1; count6 <= numaccessatom; count6++)
                //    {
                //      if (count4 == accessatom[count6])
                //      {
                //        count6 = numaccessatom + 10;
                //      }
                //    }
                //    if (count6<numaccessatom + 5)
                //    {
                //      numaccessatom++;
                //      accessatom[numaccessatom] = count4;
                //    }
                //  }
                //}
                count4 = atomnumber + 1;
                //printf(" ms 1d colapsed %4d %4d %4d icell %d framework %f distance\n", count1, count2, count3, count4, disttmp); 
              }
              else
              {
                icell_1pr[count1][count2][count3] = icell_1pr[count1][count2][count3] + tol_dist_Ar_solid/disttmp;
                if (disttmp < icelldframeatom[count1][count2][count3])
                {
                  icelldframeatom[count1][count2][count3] = disttmp;
                  icellframeatom[count1][count2][count3]  = count4;
                }
              }
            } 
          }
        }
      }
    }

    // commented in the present version that accessatoms is not in use
    //for (count1 = 1; count1 <= axesi[1]; count1++)
    //{
    //  for (count2 = 1; count2 <= axesi[2]; count2++)
    //  {
    //    for (count3 = 1; count3 <= axesi[3]; count3++)
    //    {
    //      if (icellstat[count1][count2][count3] == -3)
    //      {
    //        count4 = icellframeatom[count1][count2][count3];
    //        if (count4 != accessatom[numaccessatom])
    //        {
    //          for (count6 = 1; count6 <= numaccessatom; count6++)
    //          {
    //            if (count4 == accessatom[count6])
    //            {
    //              count6 = numaccessatom + 10;
    //            }
    //          }
    //          if (count6<numaccessatom + 5)
    //          {
    //            numaccessatom++;
    //            accessatom[numaccessatom] = count4;
    //          }  
    //        }
    //        printf(" ms 1e %4d %4d %4d icell %d nearest atom at %f distance\n", count1, count2, count3, icellframeatom[count1][count2][count3], icelldframeatom[count1][count2][count3]); 
    //      }
    //    }
    //  }
    //}

    // OJO only for test
    //for (count1 = 1; count1 <= numaccessatom; count1++)
    //{
    //  printf(" ms 1f %d accessible atom is %d\n", count1, accessatom[count1]);
    //}

    nsorticellstat = -1;
    // adding the rest of accessible atoms (dist > tol_dist_Ar_solid and also
    // determininig icellstat for non overlaping icell
    for (count1 = 1; count1 <= axesi[1]; count1++)
    {
      for (count2 = 1; count2 <= axesi[2]; count2++)
      {
        for (count3 = 1; count3 <= axesi[3]; count3++)
        {
          if (icellstat[count1][count2][count3] == -3)
          {
            if (icelldframeatom[count1][count2][count3] >= tol_dist_Ar_solid + 1.0)
            {
              icellstat[count1][count2][count3] = 26;
            }
            else
            {
              // original coding to calculate icellneighbcount
              // icellstat[count1][count2][count3] = icellneighbcount(count1, count2, count3);
              icellstat[count1][count2][count3] = 26;
            }
          }
          printf(" ms 1g %d icellstat of icell %d %d %d\n", icellstat[count1][count2][count3], count1, count2, count3);
          // original coding to consider only those close to the solid surface
          //if ( (icellstat[count1][count2][count3] < 17) && (icellstat[count1][count2][count3] >= 0) )
          //{
          //  nsorticellstat++;
          //  tsorticellstat[nsorticellstat][0] = icellstat[count1][count2][count3];
          //  tsorticellstat[nsorticellstat][1] = round(icelldframeatom[count1][count2][count3]*10000.0);
          //  tsorticellstat[nsorticellstat][2] = count1;
          //  tsorticellstat[nsorticellstat][3] = count2;
          //  tsorticellstat[nsorticellstat][4] = count3;
          //}
        }
      }
    }
    
    // original coding to consider only those close to the solid surface
    //nsorticellstat++; 
    //sorticellstat = malloc(nsorticellstat * sizeof(int*));
    //for (count1 = 0; count1 < nsorticellstat; count1++)
    //{
    //  sorticellstat[count1]    = malloc(5 * sizeof(int));
    //  sorticellstat[count1][0] = tsorticellstat[count1][0];
    //  sorticellstat[count1][1] = tsorticellstat[count1][1];
    //  sorticellstat[count1][2] = tsorticellstat[count1][2];
    //  sorticellstat[count1][3] = tsorticellstat[count1][3];
    //  sorticellstat[count1][4] = tsorticellstat[count1][4];
    //} 
    //
    //qsort(sorticellstat, nsorticellstat, sizeof sorticellstat[0], compare_qsort);

    // OJO only for test
    //for (count1 = 0; count1 < nsorticellstat; count1++)
    //{
    //  x_tmp = sorticellstat[count1][2]*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
    //  y_tmp = sorticellstat[count1][3]*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
    //  z_tmp = sorticellstat[count1][4]*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
    //  printf(" ms 1h sort_icellstat %3d %6d idist(x10000) icell %3d %3d %3d coord %.5f %.5f %.5f\n", sorticellstat[count1][0], sorticellstat[count1][1], sorticellstat[count1][2], sorticellstat[count1][3], sorticellstat[count1][4], x_tmp, y_tmp, z_tmp);
    //}

    // original coding
    // ---->>
    // put Ar atoms in all possible icells belonging to sorticellstat
    //Ar1atomnumber = 1;
    //Ar1x[Ar1atomnumber] = sorticellstat[0][2]*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
    //Ar1y[Ar1atomnumber] = sorticellstat[0][3]*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
    //Ar1z[Ar1atomnumber] = sorticellstat[0][4]*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
    //disttmp = 0.0;
    //printf(" ms 1i %d running_nsorticellstat - Ar1 insertion shortest dist %.5f icell %4d %4d %4d coord %.5f %.5f %.5f\n", Ar1atomnumber, disttmp, sorticellstat[0][2], sorticellstat[0][3], sorticellstat[0][4], Ar1x[Ar1atomnumber], Ar1y[Ar1atomnumber], Ar1z[Ar1atomnumber]);
    //for (count1 = 1; count1 < nsorticellstat; count1++)
    //{
    //  x_tmp = sorticellstat[count1][2]*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
    //  y_tmp = sorticellstat[count1][3]*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
    //  z_tmp = sorticellstat[count1][4]*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
    //  for (count2 = 1; count2 <= Ar1atomnumber; count2++)
    //  {
    //    disttmp = atomsdist(Ar1x[count2], Ar1y[count2], Ar1z[count2], x_tmp, y_tmp, z_tmp);
    //    if (disttmp < Arcell0[1]*sqrt(2.0)/2.0)
    //    {
    //      count2 = Ar1atomnumber + 10;
    //    }
    //  }
    //  // OJO only for test
    //  printf(" ms 1i1 %d running_nsorticellstat - Ar1 insertion shortest dist %.5f icell %4d %4d %4d coord %.5f %.5f %.5f\n", count1 + 1, disttmp, sorticellstat[count1][2], sorticellstat[count1][3], sorticellstat[count1][4], x_tmp, y_tmp, z_tmp);
    //  //xxxif (count2 <  Ar1atomnumber + 5)
    //  //{
    //    Ar1atomnumber++;
    //    Ar1x[Ar1atomnumber] = x_tmp; Ar1y[Ar1atomnumber] = y_tmp; Ar1z[Ar1atomnumber] = z_tmp;
    //    printf(" ms 1i2 %d running_nsorticellstat - Ar1 insertion shortest dist %.5f icell %4d %4d %4d coord %.5f %.5f %.5f\n", count1 + 1, disttmp, sorticellstat[count1][2], sorticellstat[count1][3], sorticellstat[count1][4], Ar1x[Ar1atomnumber], Ar1y[Ar1atomnumber], Ar1z[Ar1atomnumber]);
    //  //} 
    //}
    //
    ////xxx hacer prueba de poner en un fichero externo y leer luego a ver que sale
    //// OJO only for test
    //for (count1 = 1; count1 <= Ar1atomnumber; count1++)
    //{
    //  printf(" ms 1j %d Ar1 (in low space voids) %.5f %.5f %.5f\n", count1, Ar1x[count1], Ar1y[count1], Ar1z[count1]);
    //}



//xxx
    nsorticellstat = -1;
    // adding the rest of accessible atoms (dist > tol_dist_Ar_solid and also
    // determininig icellstat for non overlaping icell
    for (count1 = 1; count1 <= axesi[1]; count1++)
    {
      for (count2 = 1; count2 <= axesi[2]; count2++)
      {
        for (count3 = 1; count3 <= axesi[3]; count3++)
        {
          if (icellstat[count1][count2][count3] == -3)
          {
            if (icelldframeatom[count1][count2][count3] >= tol_dist_Ar_solid + 1.0)
            {
              icellstat[count1][count2][count3] = 26;
            }
            else
            {
              // original coding to calculate icellneighbcount
              // icellstat[count1][count2][count3] = icellneighbcount(count1, count2, count3);
              icellstat[count1][count2][count3] = 26;
            }
          }
          printf(" ms 1g %d icellstat of icell %d %d %d\n", icellstat[count1][count2][count3], count1, count2, count3);
          // original coding to consider only those close to the solid surface
          //if ( (icellstat[count1][count2][count3] < 17) && (icellstat[count1][count2][count3] >= 0) )
          //{
          //  nsorticellstat++;
          //  tsorticellstat[nsorticellstat][0] = icellstat[count1][count2][count3];
          //  tsorticellstat[nsorticellstat][1] = round(icelldframeatom[count1][count2][count3]*10000.0);
          //  tsorticellstat[nsorticellstat][2] = count1;
          //  tsorticellstat[nsorticellstat][3] = count2;
          //  tsorticellstat[nsorticellstat][4] = count3;
          //}
        }
      }
    }

    
    count6 = 0; count7 = 0; 
    tol_dist_Ar_Ar = 0.1*dist_Ar_Ar;
    for (count1 = 1; count1 <= Ar1atomnumber; count1++) // new coding for deal with Ar1
    {
      //Arxca[count1] = Ar1x[count1]; Aryca[count1] = Ar1y[count1]; Arzca[count1] = Ar1z[count1];
      removeAr[count1] = 0;
      printf(" ms 12c1 %4d running_Ar1atomnumber - coord %.5f %.5f %.5f Ar1 deleted? %4d removed Ar %4d from %4d\n", count1,  Ar1x[count1], Ar1y[count1], Ar1z[count1], removeAr[count1], count6, Ar1atomnumber);
      disttmp1 = 10000.0; // OJO eliminar on standard running
      for (count2 = 1; count2 <= Aratomnumber; count2++)
      {
        disttmp = atomsdist(Ar1x[count1], Ar1y[count1], Ar1z[count1], Arx[count2], Ary[count2], Arz[count2]);
        if (disttmp < disttmp1)
        {
          disttmp1 = disttmp;
        }
        if (disttmp < tol_dist_Ar_Ar) // aqui cambiar el factor
        {
          removeAr[count1] = 1;
          count2 = Aratomnumber + 10;
          count6++;
        }
      }
      // OJO only for test
      printf(" ms 12c2 %4d running_Ar1atomnumber - coord %.5f %.5f %.5f Ar1 deleted? %4d removed Ar %4d from %4d dist_min %.5f\n", count1,  Ar1x[count1], Ar1y[count1], Ar1z[count1], removeAr[count1], count6, Ar1atomnumber, disttmp1);
      if (removeAr[count1] == 0)
      {
        for (count2 = 1;  count2 < count1; count2++)
        {
          if (removeAr[count2] == 0)
          {
            disttmp = atomsdist(Ar1x[count1], Ar1y[count1], Ar1z[count1], Ar1x[count2], Ar1y[count2], Ar1z[count2]);
            if (disttmp < Arcell0[1]*sqrt(2.0)/2.0)
            {
              removeAr[count1] = 1;
              count2 = count1 + 10;
              count7++;
            }
          }
        }
        // OJO only for test
        printf(" ms 12d %d running_Ar1atomnumber 2 - Ar1 deleted? %d insertion shortest dist %.5f removed Ar %d from %d\n", count1, removeAr[count1], disttmp, count7, Ar1atomnumber);
      }
    }
   
    // only for checking debuging purpose, comment on standard running
    count4 = maxloading;
    for (count1 = 1; count1 <= Ar1atomnumber; count1++)
    {
      if (removeAr[count1] == 0)
      {
        count4++;
        printf(" ms12e Ar %d\t Loading Min-Max %d-%d Coord cryst %.5f %.5f %.5f Ar1 number %d\n", count4, minloading, count4, Ar1x[count1], Ar1y[count1], Ar1z[count1], Ar1atomnumber);
      }
    }
    printf("\n");
// xxx aqui
    for (count1 = 1; count1 <= Ar1atomnumber; count1++)
    {
      if (removeAr[count1] == 0)
      {
        count3++;
        fprintf(OutpF, " O%d O %.5f %.5f %.5f 0.00\n", count3, Ar1x[count1], Ar1y[count1], Ar1z[count1]);
      }
    }
*/
    fclose(OutpF);
    
}   /***********    end of MainJob     ************/

/****************************************************************************
     This function counts the number of available cells in the direct 
     neighbourhood of the current icell. Note that there are 26 neighbours
****************************************************************************/
int    icellneighbcount(int ix, int iy, int iz)
{
    int    count1, count2, count3, count4, licellneighbcount;
    int    x_m1, x_p1, y_m1, y_p1, z_m1, z_p1;    

    licellneighbcount = 0; 
    x_m1 = ix -1; x_p1 = ix + 1;  
    y_m1 = iy -1; y_p1 = iy + 1;
    z_m1 = iz -1; z_p1 = iz + 1;

    if (ix == 1)
    {
      x_m1 = axesi[1];
    }
    else
    {
      if (ix == axesi[1])
      {
        x_p1 = 1;
      }
    }
    if (iy == 1)
    {
      y_m1 = axesi[2];
    }
    else
    {
      if (iy == axesi[2])
      {
        y_p1 = 1;
      }
    }
    if (iz == 1)
    {
      z_m1 = axesi[3];
    }
    else
    {
      if (iz == axesi[3])
      {
        z_p1 = 1;
      }
    }

    for (count4 = 1; count4 <= 26; count4++)
    {
      switch (count4)
      {
        case 1:
          count1 = x_m1; count2 = y_p1; count3 = z_m1;
        break;
        case 2:
          count1 = ix;   count2 = y_p1; count3 = z_m1;
        break;
        case 3:
          count1 = x_p1; count2 = y_p1; count3 = z_m1;
        break;
        case 4:
          count1 = x_m1; count2 = iy;   count3 = z_m1;
        break;
        case 5:
          count1 = ix;   count2 = iy;   count3 = z_m1;
        break;
        case 6:
          count1 = x_p1; count2 = iy;   count3 = z_m1;
        break;
        case 7:
          count1 = x_m1; count2 = y_m1; count3 = z_m1;
        break;
        case 8:
          count1 = ix;   count2 = y_m1; count3 = z_m1;
        break;
        case 9:
          count1 = x_p1; count2 = y_p1; count3 = z_m1;
        break;
	case 10:
          count1 = x_m1; count2 = y_p1; count3 = iz;  
        break;
	case 11:
          count1 = ix;   count2 = y_p1; count3 = iz;
        break;
	case 12:
          count1 = x_p1; count2 = y_p1; count3 = iz;
        break;
        case 13:
          count1 = x_m1; count2 = iy;   count3 = iz;
        break;
        case 14:
          count1 = x_p1; count2 = iy;   count3 = iz;
        break;
        case 15:
          count1 = x_m1; count2 = y_m1; count3 = iz;
        break;
        case 16:
          count1 = ix;   count2 = y_m1; count3 = iz;
        break;
        case 17:
          count1 = x_p1; count2 = y_m1; count3 = iz;
        break;
	case 18:
          count1 = x_m1; count2 = y_p1; count3 = z_p1;
        break;
	case 19:
          count1 = ix;   count2 = y_p1; count3 = z_p1;
        break;
	case 20:
          count1 = x_p1; count2 = y_p1; count3 = z_p1;
        break;
        case 21:
          count1 = x_m1; count2 = iy;   count3 = z_p1;
        break;
        case 22:
          count1 = ix;   count2 = iy;   count3 = z_p1;
        break;
        case 23:
          count1 = x_p1; count2 = iy;   count3 = z_p1;
        break;
        case 24:
          count1 = x_m1; count2 = y_m1; count3 = z_p1;
        break;
        case 25:
          count1 = ix;   count2 = y_m1; count3 = z_p1;
        break;
        case 26:
          count1 = x_p1; count2 = y_m1; count3 = z_p1;
        break;
      } 
      
      if ( (icellstat[count1][count2][count3] == -3) || (icellstat[count1][count2][count3] >= 0) )
      {
        licellneighbcount++;
      }  
    }

    return licellneighbcount;

}   /***********   end of icellneighbcount ************/    


/**************************************************************************
qsort function tested by Alfre from internet
**************************************************************************/
int    compare_qsort( const void *pa, const void *pb)
//int    compare_qsort( const void *pa, const void *pb, const void *pc, const void *pd)
{
    const int *a = *(const int **)pa;
    const int *b = *(const int **)pb;
    if(a[0] == b[0])
        return a[1] - b[1];
    else
        return a[0] - b[0];
}   /***********   end of compare_qsort   ************/


/**************************************************************************
  label icellstat as -2 or -1 by overlaping with framework atoms or Ar,
  respectiveli
  Ar_or_atom = 1 framework and 2 Ar
**************************************************************************/
void   icelloverlap_Ar_atom(int Ndisp_a1, int Ndisp_b1, int Ndisp_c1, int Ar_or_atom, double xatom1, double yatom1, double zatom1)
{
    int    count1, count2, count3, count4, new_icellstat;
    int    count5, count6, count7, count8, x_i, y_i, z_i;
    double tmp_tol, x11, y11, z11, disttmp;
    int    disp_count1, disp_count2, disp_count3;

    // x_i = xatom1*axesidelta[1];
    // y_i = yatom1*axesidelta[2];
    // z_i = zatom1*axesidelta[3];

    x_i = xatom1*axesi[1] + 1;  // icell_x where the running atom belongs
    y_i = yatom1*axesi[2] + 1;
    z_i = zatom1*axesi[3] + 1;

    if (Ar_or_atom == 1)
    {
      new_icellstat = -2;
      tmp_tol       = tol_dist_Ar_solid;
    }
    else
    { 
      if (Ar_or_atom == 2)
      {
        new_icellstat = -1;
        tmp_tol       = tol_dist_Ar_Ar;
      }
      else
      {
        new_icellstat = -4;
        tmp_tol       = tol_dist_Ar_Ar;
      }
    }

    //printf(" ms icelloverlap_Ar_atom_1 type %d icell %3d  %3d %3d label_icellstat %d tolerance %.5f\n", Ar_or_atom, x_i, y_i, z_i, new_icellstat, tmp_tol);
 
   // for (disp_count1 = x_i - Ndisp_a1; disp_count1 <= x_i + Ndisp_a1; disp_count1++)
   // {
   //   x11 = x_i*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
   //   for (disp_count2 = y_i - Ndisp_b1; disp_count2 <= y_i + Ndisp_b1; disp_count2++)
   //   { 
   //     y11 = y_i*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
   //     for (disp_count3 = z_i - Ndisp_c1; disp_count3 <= z_i + Ndisp_c1; disp_count3++)
   //     {
   //       if (icellstat[count5][count6][count7] < -2)
   //       {   
   //         z11 = z_i*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
   //         disttmp = atomsdist(xatom1, yatom1, zatom1, x11, y11, z11);
   //         if (disttmp < tmp_tol)
   //         {
   //           icellstat[count5][count6][count7] = new_icellstat;
   //         }
   //         printf(" ms icelloverlap_Ar_atom_2 type %d %.5f icell_double_dist lower_double_tol %.5f\n", Ar_or_atom, disttmp, tmp_tol);
   //       }
   //     }  
   //   }    
   // }


   for (disp_count1 = x_i - Ndisp_a1; disp_count1 <= x_i + Ndisp_a1; disp_count1++)
   {
     if (disp_count1 <= 0)
     {
       count5 = axesi[1] + disp_count1 + 1;
     }
     else
     {
       if (disp_count1 > axesi[1])
       {
         count5 = disp_count1 - axesi[1];
       }
       else
       {
         count5 = disp_count1;
       }
     }
     x11 = count5*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
     //count1 = (x_i - disp_count1)*(x_i - disp_count1);
     for (disp_count2 = y_i - Ndisp_b1; disp_count2 <= y_i + Ndisp_b1; disp_count2++)
     { 
       if (disp_count2 <= 0)
       {
         count6 = axesi[2] + disp_count2 + 1;
       }
       else
       {
         if (disp_count2 > axesi[2])
         {
           count6 = disp_count2 - axesi[2];
         }
         else
         {
           count6 = disp_count2;
         }
       }
       y11 = count6*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
       //count2 = (y_i - disp_count2)*(y_i - disp_count2);
       for (disp_count3 = z_i - Ndisp_c1; disp_count3 <= z_i + Ndisp_c1; disp_count3++)
       {
         if (disp_count3 <= 0)
         {
           count7 = axesi[3] + disp_count3 + 1;
         }
         else
         {
           if (disp_count3 > axesi[3])
           {
             count7 = disp_count3 - axesi[3];
           }
           else
           {
             count7 = disp_count3;
           }
         } 
         z11 = count7*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
         //count3 = (z_i - disp_count3)*(z_i - disp_count3);
        
         if (icellstat[count5][count6][count7] < -2)
         {   
           disttmp = atomsdist(xatom1, yatom1, zatom1, x11, y11, z11);
           if (disttmp < tmp_tol)
           {
             icellstat[count5][count6][count7] = new_icellstat;
           }
           //printf(" ms icelloverlap_Ar_atom_2 type %d icell  %3d %3d %3d  %.5f icell_double_dist lower_double_tol %.5f  icellstat %d\n", Ar_or_atom, count5, count6, count7, disttmp, tmp_tol, icellstat[count5][count6][count7]);
         }
       }  
     }    
   }




   // original coding using counts of icells in some cases instead of calculated distances
   //
      
   //for (disp_count1 = x_i - Ndisp_a1; disp_count1 <= x_i + Ndisp_a1; disp_count1++)
   //{
   //  if (disp_count1 <= 0)
   //  {
   //    count5 = axesi[1] + disp_count1 + 1;
   //  }
   //  else
   //  {
   //    if (disp_count1 > axesi[1])
   //    {
   //      count5 = disp_count1 - axesi[1];
   //    }
   //    else
   //    {
   //      count5 = disp_count1;
   //    }
   //  }
   //  count1 = (x_i - disp_count1)*(x_i - disp_count1);
   //  for (disp_count2 = y_i - Ndisp_b1; disp_count2 <= y_i + Ndisp_b1; disp_count2++)
   //  { 
   //    if (disp_count2 <= 0)
   //    {
   //      count6 = axesi[2] + disp_count2 + 1;
   //    }
   //    else
   //    {
   //      if (disp_count2 > axesi[2])
   //      {
   //        count6 = disp_count2 - axesi[2];
   //      }
   //      else
   //      {
   //        count6 = disp_count2;
   //      }
   //    }
   //    count2 = (y_i - disp_count2)*(y_i - disp_count2);
   //    for (disp_count3 = z_i - Ndisp_c1; disp_count3 <= z_i + Ndisp_c1; disp_count3++)
   //    {
   //      if (disp_count3 <= 0)
   //      {
   //        count7 = axesi[3] + disp_count3 + 1;
   //      }
   //      else
   //      {
   //        if (disp_count3 > axesi[3])
   //        {
   //          count7 = disp_count3 - axesi[3];
   //        }
   //        else
   //        {
   //          count7 = disp_count3;
   //        }
   //      } 
   //      count3 = (z_i - disp_count3)*(z_i - disp_count3);
   //     
   //      if (icellstat[count5][count6][count7] < -2)
   //      {   
   //        count8 = round(tmp_tol);
   //        count8 = count8 * count8;
   //        if ( (count1 + count2 + count3) <= count8 )
   //        {
   //          icellstat[count5][count6][count7] = new_icellstat;
   //          printf(" ms icelloverlap_Ar_atom_2a type %d %3d icell_int_dist^2 lower_int_tol %d\n", Ar_or_atom, count1 + count2 + count3, count8); 
   //        }   
   //        else
   //        {
   //          x11 = x_i*1.0/(axesi[1]*1.0) + 1.0/(2.0*axesi[1]);
   //          y11 = y_i*1.0/(axesi[2]*1.0) + 1.0/(2.0*axesi[2]);
   //          z11 = z_i*1.0/(axesi[3]*1.0) + 1.0/(2.0*axesi[3]);
   //          disttmp = atomsdist(xatom1, yatom1, zatom1, x11, y11, z11);
   //          if (disttmp < tmp_tol)
   //          {
   //            icellstat[count5][count6][count7] = new_icellstat;
   //          }
   //          printf(" ms icelloverlap_Ar_atom_2b type %d %.5f icell_double_dist lower_double_tol %.5f\n", Ar_or_atom, disttmp, tmp_tol);
   //        }
   //      }
   //    }  
   //  }    
   //}
}   /***** ******   end of icelloverlap_Ar_atom   ************/
           
           
/********* *******************************************************************
     this  function extractes the atoms names, core-shel specification
     and c oordinates
           
     types of files allowed
     1:   gulp input and/or restar   *.gin, *.res
     2:   gulp output                *.gout or *.gt
     3:   msi xtl                    *.xtl
     4:   gulp car file              *.car  NOTE, is GULP CAR FILE
     5:   cif                        *.cif

     the results are storaged in external variables:
     double xatom[MaxAtomNum], yatom[MaxAtomNum], zatom[MaxAtomNum];
     char   atomname[MaxAtomNum][4];
     char   atom_c_s[MaxAtomNum];
     int    atomnumber;
     double cell[N];
*****************************************************************************/
int    getatomcoor(char coordFName[Maxline], int CoorFType)
{
        double r1, r2, r3, r4, tempval1;
        int    c1, c2, c3, c4, tmp_dummy, cont1, maxline, intspace, contspace;
        char   s1[Maxline], s2[Maxline], s3[Maxline], s4[Maxline];
        char   s5[Maxline], s6[Maxline], s7[Maxline], s8[Maxline];
        char   line[Maxline];
        int    readinline;
        int    getatomcoor;
        FILE  *InProcesFile;
        int    count1, count2, count3;
        char  *chtmp1;
        char   s1a[Maxline];

        //printf("already in getatomcoor 0, reading file %s of type %d \n", coordFName, CoorFType);
        
        InProcesFile = fopen(coordFName,"r");
        getatomcoor = 1; 
        maxline    = Maxline;
        atomnumber = 0;
        readinline = 1;
        
        //printf("already in getatomcoor 1, reading file %s of type %d \n", coordFName, CoorFType);
        
        switch (CoorFType)
        {
        case 1:
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "cell", 4)!= 0)
           {
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           fgets(line, maxline, InProcesFile);
           /* blanck line skipper */
           if (line == NULL)
             {
             getatomcoor = 0;
             break;
             }
           while (blanckline(line) == 0)
              {
                fgets(line, maxline, InProcesFile);
              }
           /*  end of blanck line skipper  */
           sscanf(line,"%lf %lf %lf %lf %lf %lf", &cell[1], &cell[2], &cell[3], &cell[4], &cell[5], &cell[6]);
           strcpy(s1, "fractional");
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
           while (readinline == 1)
             {
               fgets(line, maxline, InProcesFile);
               if (line == NULL)
                {
                  getatomcoor = 0;
                  break;
                }
               while (blanckline(line) == 0)
                {
                  fgets(line, maxline, InProcesFile);
                }

/*               printf(" aviso0 %s\n", line);      */
               if (line != NULL)
                 {
                   if ( countwords(line, 0,strlen(line)-1) <= 3)
                     {
                       readinline = 0;
                     }
                     else
                     {
/*                       printf(" aviso1 %d\n", readinline);     */
                       tmp_dummy = 1;
                       strcpy(s2 , grepnword(line, 0, 1));
/*                       printf(" palabra1 %s\n", s2);           */
                       strcpy(s3 , grepnword(line, 0, 2));
/*                       printf(" palabra2 %s\n", s3);           */
                       strcpy(s4 , grepnword(line, 0, 3));
/*                       printf(" palabra3 %s\n", s4);           */
                       strcpy(s5 , grepnword(line, 0, 4));
/*                       printf(" palabra4 %s\n", s5);           */
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            if ( countwords(line, 0,strlen(line)-1) <= 4)
                              readinline = 0;
                            else
                              strcpy(s6 , grepnword(line, 0, 5));
/*                            printf(" palabra5 %s\n", s6);      */
                            if (readinline != 0)
                               {
/*                                 printf(" aviso s6a %d\n", readinline);   */
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     {
                                       tmp_dummy = 1;
                                     }
                                   else
                                     {
                                       readinline = 0;
                                     }
                               }
                               else
                                 {
                                   readinline = 0;
/*                                   printf(" aviso s6b\n");      */
                                 }
                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {  
                                 readinline = 0;
                               }
                               
                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {
                                 readinline = 0;
                               }
                                                          
/*                            printf(" aviso readinline = %d\n", readinline);       */
                            if (readinline != 0)
                              {
                                atomnumber++;
                                sscanf(s4, "%lf", &xatom[atomnumber-1]);
/*                                printf("%f\t", xatom[atomnumber-1]);   */
                                sscanf(s5, "%lf", &yatom[atomnumber-1]);
/*                                printf(" %f\t", yatom[atomnumber-1]);  */
                                sscanf(s6, "%lf", &zatom[atomnumber-1]);
/*                                printf("%f\n", zatom[atomnumber-1]);   */
                                
                                tmp_dummy = 0;
                                for (tmp_dummy = 0; tmp_dummy <= strlen(s2) -1; tmp_dummy++)
                                   atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy];                                
                                /*
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                */
                                tmp_dummy = 1;
/*                                printf(" nombre atomo ?? %d %s\n", atomnumber-1, atomname[atomnumber-1]);  */
                               
                                atom_c_s[atomnumber-1] = s3[0];
                              }
                          }
                          else
                          {
                            printf(" llega a aqui ???");            
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
/*                          printf(" retorno de la columna 27\n");  */
                          }
/*                     printf(" retorno de la columna 22\n");       */
                     }
/*                 printf(" retorno de la columna 18\n");           */
                 }
                 else
                   readinline = 0;  
/*             printf(" retorno de la columna 14\n");               */
             }
             
             printf("already in getatomcoor 2\n");

           break;

        case 2:
           
           strcpy(s1, "Final fractional coordinates of atoms");
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
           for (c1 = 1; c1 <= 6; c1++)
               {
               fgets(line, maxline, InProcesFile);
               }
           while (readinline == 1)
             {
               sscanf(line,"%s", s1);
               while ( (line != NULL) && (isspace(s1[0]) == 0) )
                  {
                  fgets(line, maxline, InProcesFile);
                  sscanf(line,"%s", s1);
                  }

               if (line != NULL)
                 {
                   if ( (strncmp(line, "--------------------", 20) != 0) || (strncmp(line, "Final cell", 10)!= 0) )
                     {
                       readinline = 0;
                     }
                     else
                     {
                       tmp_dummy = 1;
                       strcpy(s2 , grepnword(line, 0, 2));
                       strcpy(s3 , grepnword(line, 0, 3));
                       strcpy(s4 , grepnword(line, 0, 4));
                       strcpy(s5 , grepnword(line, 0, 5));
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            strcpy(s6 , grepnword(line, 0, 6));
                            if (strstr(s6, "") == NULL)
                               {
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     tmp_dummy = 1;
                                   else
                                     readinline = 0;
                               }
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf(s3[0], "%c", atom_c_s[atomnumber-1]);
                                sprintf(s4, "%.5f",xatom[atomnumber-1]);
                                sprintf(s5, "%.5f",yatom[atomnumber-1]);
                                sprintf(s6, "%.5f",zatom[atomnumber-1]);
                              }
                          }
                          else
                          {
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
                          }
                     }
                 }
                 else
                   readinline = 0;  
             }
           
           strcpy(s1, "Final cell parameters and derivatives");
           tmp_dummy = 1000;
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
                                      
           if (tmp_dummy != 0)
           {
             for (c1 = 1; c1 <= 6; c1++)
                 cell[c1] = -1.0;
                 c3 = 1; 
           }                                                                                                                                       
           else
           {
             for (c1 = 1; c1 <= 3; c1++)
               {
               fgets(line, maxline, InProcesFile);
               }
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[1], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s", s1, &cell[2], s2);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[3], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[4], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[5], s2, s3, &tempval1, s4);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s %lf %s %s %lf %s", s1, &cell[6], s2, s3, &tempval1, s4);
           }
           break;

        case 3:
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "CELL", 4)!= 0)
           {   
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           fgets(line, maxline, InProcesFile);
           /* blanck line skipper */
           if (line == NULL)
             {
             getatomcoor = 0;
             break;
             }
           while (blanckline(line) == 0)
              {
                fgets(line, maxline, InProcesFile);
              }
           /*  end of blanck line skipper  */
           sscanf(line,"%lf %lf %lf %lf %lf %lf", &cell[1], &cell[2], &cell[3], &cell[4], &cell[5], &cell[6]);
           strcpy(s1, "NAME");
           tmp_dummy = strgrepf(line, s1, maxline, InProcesFile);
/*           fgets(line, maxline, InProcesFile);   */
           while (readinline == 1)
             {
               fgets(line, maxline, InProcesFile);
               if (line == NULL)
                {
                  getatomcoor = 0;
                  break;
                }
               while (blanckline(line) == 0)
                {
                  fgets(line, maxline, InProcesFile);
                }

/*               printf(" aviso0 %s\n", line);     */
               if (line != NULL)
                 {
                   if ( countwords(line, 0,strlen(line)-1) <= 3)
                     {
                       readinline = 0;
                     }
                     else
                     {
/*                       printf(" aviso1 %d\n", readinline);     
                       tmp_dummy = 1;
                       strcpy(s2 , grepnword(line, 0, 1));
                      printf(" palabra1 %s\n", s2);           
                       strcpy(s3 , "c");
                       printf(" palabra2 %s\n", s3);           
                       strcpy(s4 , grepnword(line, 0, 2));
                       printf(" palabra3 %s\n", s4);           
                       strcpy(s5 , grepnword(line, 0, 3));
*/
                       sscanf(line,"%s %s %s %s", s2, s4, s5, s6);
                       strcpy(s3 , "c");
/*                       printf(" read %s %s %s %s\n", s2, s3, s4, s5);          */
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            if ( countwords(line, 0,strlen(line)-1) <= 3)
                              readinline = 0;
                            else
                               tmp_dummy = 1;
/*                              strcpy(s6 , grepnword(line, 0, 4)); */
/*                            printf(" palabra5 %s\n", s6);      */
                            if (readinline != 0)
                               {
/*                                 printf(" aviso s6a %d\n", readinline);   */
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     {
                                       tmp_dummy = 1;
                                     }
                                   else
                                     {
                                       readinline = 0;
                                     }
                               }
                               else
                                 {
                                   readinline = 0;
/*                                   printf(" aviso s6b\n");      */
                                 }
                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {  
                                 readinline = 0;
                               }
                               
                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {
                                 readinline = 0;
                               }
                                                          
/*                            printf(" aviso readinline = %d\n", readinline);       */
                            if (readinline != 0)
                              {
                                atomnumber++;
                                sscanf(s4, "%lf", &xatom[atomnumber-1]);
/*                                printf("%f\t", xatom[atomnumber-1]);   */
                                sscanf(s5, "%lf", &yatom[atomnumber-1]);
/*                                printf(" %f\t", yatom[atomnumber-1]);  */
                                sscanf(s6, "%lf", &zatom[atomnumber-1]);
/*                                printf("%f\n", zatom[atomnumber-1]);   */
                                
                                tmp_dummy = 0;
                                for (tmp_dummy = 0; tmp_dummy <= strlen(s2) -1; tmp_dummy++)
                                   atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy];                                
                                /*
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                */
                                tmp_dummy = 1;
/*                                printf(" nombre atomo ?? %d %s\n", atomnumber-1, atomname[atomnumber-1]);  */
                               
                                atom_c_s[atomnumber-1] = s3[0];
/*                                printf(" %c\n", atom_c_s[atomnumber-1]);                               
                                                              
                                printf("ojo %d\t %s\t %c\t", atomnumber-1,  atomname[atomnumber-1], atom_c_s[atomnumber-1]);
                                printf("%f\t %f\t %f; %s\n", xatom[atomnumber-1], yatom[atomnumber-1], zatom[atomnumber-1], line); 
*/
                              }
                          }
                          else
                          {
                            printf(" llega a aqui ???");            
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
/*                          printf(" retorno de la columna 27\n");  */
                          }
/*                     printf(" retorno de la columna 22\n");       */
                     }
/*                 printf(" retorno de la columna 18\n");           */
                 }
                 else
                   readinline = 0;  
/*             printf(" retorno de la columna 14\n");               */
             }
             
             printf("already in getatomcoor 2\n");

           break;

        case 4:
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "PBC", 3)!= 0)
           {   
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           while (strncmp(s1, "PBC", 3)!= 0)
           {
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }

           if (line == NULL)
             {
             getatomcoor = 0;
             break;
             }

/* need to put a counter, to stop while loop after reaching 100 000 lines in both readings above */
 
           sscanf(line,"%s %lf %lf %lf %lf %lf %lf", s1, &cell[1], &cell[2], &cell[3], &cell[4], &cell[5], &cell[6]);

           /* blanck line skipper */
           while (blanckline(line) == 0)
              {
                fgets(line, maxline, InProcesFile);
              }
           /*  end of blanck line skipper  */

           while (readinline == 1)
             {
               fgets(line, maxline, InProcesFile);
               if (line == NULL)
                {
                  getatomcoor = 0;
                  break;
                }
               while (blanckline(line) == 0)
                {
                  fgets(line, maxline, InProcesFile);
                }

/*               printf(" aviso0 %s\n", line);     */
               if (line != NULL)
                 {
                   if ( countwords(line, 0,strlen(line)-1) <= 3)
                     {
                       readinline = 0;
                     }
                     else
                     {
/*                       printf(" aviso1 %d\n", readinline);     
                       tmp_dummy = 1;
                       strcpy(s2 , grepnword(line, 0, 1));
                      printf(" palabra1 %s\n", s2);           
                       strcpy(s3 , "c");
                       printf(" palabra2 %s\n", s3);           
                       strcpy(s4 , grepnword(line, 0, 2));
                       printf(" palabra3 %s\n", s4);           
                       strcpy(s5 , grepnword(line, 0, 3));
*/
                       sscanf(line,"%s %s %s %s", s2, s4, s5, s6);
                       strcpy(s3 , "c");
/*                       printf(" read %s %s %s %s\n", s2, s3, s4, s5);          */
                       if ( (s3[0] == 'c') || (s3[0] == 's') )
                          {
                            if ( countwords(line, 0,strlen(line)-1) <= 3)
                              readinline = 0;
                            else
                               tmp_dummy = 1;
/*                              strcpy(s6 , grepnword(line, 0, 4)); */
/*                            printf(" palabra5 %s\n", s6);      */
                            if (readinline != 0)
                               {
/*                                 printf(" aviso s6a %d\n", readinline);   */
                                 if ( (s6[0] =='+') || (s6[0] =='-') || (s6[0] =='.') || (isdigit(s6[0]) != 0) )
                                     {
                                       tmp_dummy = 1;
                                     }
                                   else
                                     {
                                       readinline = 0;
                                     }
                               }
                               else
                                 {
                                   readinline = 0;
/*                                   printf(" aviso s6b\n");      */
                                 }
                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {  
                                 readinline = 0;
                               }
                               
                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                               {
                                 tmp_dummy = 1;
                               }
                               else
                               {
                                 readinline = 0;
                               }
                                                          
/*                            printf(" aviso readinline = %d\n", readinline);       */
                            if (readinline != 0)
                              {
                                atomnumber++;
                                sscanf(s4, "%lf", &xatom[atomnumber-1]);
/*                                printf("%f\t", xatom[atomnumber-1]);   */
                                sscanf(s5, "%lf", &yatom[atomnumber-1]);
/*                                printf(" %f\t", yatom[atomnumber-1]);  */
                                sscanf(s6, "%lf", &zatom[atomnumber-1]);
/*                                printf("%f\n", zatom[atomnumber-1]);   */
                                
                                tmp_dummy = 0;
                                for (tmp_dummy = 0; tmp_dummy <= strlen(s2) -1; tmp_dummy++)
                                   atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy];                                
                                /*
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                */
                                tmp_dummy = 1;
/*                                printf(" nombre atomo ?? %d %s\n", atomnumber-1, atomname[atomnumber-1]);  */
                               
                                atom_c_s[atomnumber-1] = s3[0];
/*                                printf(" %c\n", atom_c_s[atomnumber-1]);                               
                                                              
                                printf("ojo %d\t %s\t %c\t", atomnumber-1,  atomname[atomnumber-1], atom_c_s[atomnumber-1]);
                                printf("%f\t %f\t %f; %s\n", xatom[atomnumber-1], yatom[atomnumber-1], zatom[atomnumber-1], line); 
*/
                              }
                          }
                          else
                          {
                            printf(" llega a aqui ???");            
                            if ( (s3[0] =='+') || (s3[0] =='-') || (s3[0] =='.') || (isdigit(s3[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s4[0] =='+') || (s4[0] =='-') || (s4[0] =='.') || (isdigit(s4[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if ( (s5[0] =='+') || (s5[0] =='-') || (s5[0] =='.') || (isdigit(s5[0]) != 0) )
                                 tmp_dummy = 1;
                               else
                                 readinline = 0;

                            if (readinline != 0)
                              {
                                atomnumber++;
                                tmp_dummy = 1;
                                while ( (atomname[atomnumber-1][tmp_dummy] = s2[tmp_dummy]) != '\0')
                                  tmp_dummy++;
                                tmp_dummy = 1;
                                sprintf('c', "%c", atom_c_s[atomnumber-1]);
                                sprintf(s3, "%.5f",xatom[atomnumber-1]);
                                sprintf(s4, "%.5f",yatom[atomnumber-1]);
                                sprintf(s5, "%.5f",zatom[atomnumber-1]);
                              }
/*                          printf(" retorno de la columna 27\n");  */
                          }
/*                     printf(" retorno de la columna 22\n");       */
                     }
/*                 printf(" retorno de la columna 18\n");           */
                 }
                 else
                   readinline = 0;  
/*             printf(" retorno de la columna 14\n");               */
             }
             
             printf("already in getatomcoor 2\n");

           break;

        case 5:
           fgets(line, maxline, InProcesFile);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           sscanf(line,"%s", s1);
           while (strncmp(s1, "_cell_length_a", 14)!= 0)
           {
//             printf("%s", line);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           memset(s1, '\0', sizeof(s1));
           memset(s2, '\0', sizeof(s2));
           //printf("test reading line _cell_length_a %s", line);
           sscanf(line,"%s %s", s1, s2);
           //printf("test reading text _cell_length_a %s\n", s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3 + 10;
             }
           }
           memset(s1a, '\0', sizeof(s1a));
           if ( count1 >= (count3 + 10) )
           {
             strcpy(s1a, grepnmword(s2, 0, count2));
           }
           else
           {
             strcpy(s1a, s2);
           }
           //strcpy(s1a, grepnmword(s2, 0, count2));
           //sscanf(grepnmword(s2, 0, count2), "%s", s1a);
           // printf("test reading s2 %s s1a %s and determining count2 %d\n", s2, s1a, count2);
           sscanf(s1a, "%lf", &cell[1]);
           // temporal para este problema!!!
           // sscanf(s2, "%lf", &cell[1]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           //printf("test reading line _cell_length_b %s", line);
           sscanf(line,"%s %s", s1, s2);
           //printf("test reading text _cell_length_b %s\n", s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3 + 10;
             }
           }
           memset(s1, '\0', sizeof(s1));
           if ( count1 >= (count3 + 10) )
           {
             strcpy(s1, grepnmword(s2, 0, count2));
           }
           else
           {
             strcpy(s1, s2);
           }
           sscanf(s1, "%lf", &cell[2]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           //printf("test reading line _cell_length_c %s", line);
           sscanf(line,"%s %s", s1, s2);
           //printf("test reading text _cell_length_c %s\n", s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3 + 10;
             }
           }
           memset(s1, '\0', sizeof(s1));
           if ( count1 >= (count3 + 10) )
           {
             strcpy(s1, grepnmword(s2, 0, count2));
           }
           else
           {
             strcpy(s1, s2);
           }
           sscanf(s1, "%lf", &cell[3]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile); 
           //printf("test reading line _cell_length_alpha %s", line);
           sscanf(line,"%s %s", s1, s2);
           //printf("test reading text _cell_length_alpha %s\n", s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3 + 10;
             }
           }
           memset(s1, '\0', sizeof(s1));
           if ( count1 >= (count3 + 10) )
           {
             strcpy(s1, grepnmword(s2, 0, count2));
           }
           else
           {
             strcpy(s1, s2);
           }
           sscanf(s1, "%lf", &cell[4]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           //printf("test reading line _cell_length_beta %s", line);
           sscanf(line,"%s %s", s1, s2);
           //printf("test reading text _cell_length_beta %s\n", s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3 + 10;
             }
           }
           memset(s1, '\0', sizeof(s1));
           if ( count1 >= (count3 + 10) )
           {
             strcpy(s1, grepnmword(s2, 0, count2));
           }
           else
           {
             strcpy(s1, s2);
           }
           sscanf(s1, "%lf", &cell[5]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           //printf("test reading line _cell_length_gamma %s", line);
           sscanf(line,"%s %s", s1, s2);
           //printf("test reading text _cell_length_gamma %s\n", s2);
           count3 = strlen(s2) - 1;
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3 + 10;
             }
           }
           memset(s1, '\0', sizeof(s1));
           if ( count1 >= (count3 + 10) )
           {
             strcpy(s1, grepnmword(s2, 0, count2));
           }
           else
           {
             strcpy(s1, s2);
           }
           sscanf(s1, "%lf", &cell[6]);
           printf("Read cell parameters %f %f %f %f %f %f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);
           
           while (strncmp(s1, "_atom_site_", 11) != 0)
           {
             chtmp1 = fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           //printf("1 atom_site %s\n", s1);
           chtmp1 = fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           //printf("2 atom_site %s\n", s1);
           tmp_dummy = 1;
           if (strncmp(s1, "_atom_site_fract_x", 18)!= 0)
           {
           	 tmp_dummy = 2;
		   } 
		   //printf(" number of atom label columns %d\n", tmp_dummy);      
		   
           while (strncmp(s1, "_atom_site_fract_z", 18)!= 0)
           {
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           while (strncmp(s1, "_atom_site_", 11) == 0)
           {
             chtmp1 = fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           while ( (chtmp1 != NULL) && (strncmp(s1, "loop_", 5) != 0) )
           {
             if (strlen(line)>2)
             {
             	if (tmp_dummy == 2)
             	{
                  sscanf(line,"%s %s %lf %lf %lf %lf", s2, atomname[atomnumber], &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);
             	}
             	else
             	{
             	  sscanf(line,"%s %lf %lf %lf %lf", atomname[atomnumber], &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);  
             	}
                //printf("%d %s %.6f %.6f %.6f %.6f type input %d\n", atomnumber, atomname[atomnumber], xatom[atomnumber], yatom[atomnumber], zatom[atomnumber], atomcharge[atomnumber], tmp_dummy);
                atomnumber++;
             }
             chtmp1 = fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }

           strcpy(s3 , "c");
           for (count1 = 0; count1 < atomnumber; count1++)
           {
             atom_c_s[count1] = s3[0];
           }

           //printf("already in getatomcoor 2\n");

           break;
        }
      
         printf("atomnumber %d\n", atomnumber); 
         //printf("cell  %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]); 
         //for (c1 = 0; c1 <= atomnumber-1; c1++)
         //      {
         //        printf("%s %c %.5f %.5f %.5f\n", atomname[c1], atom_c_s[c1], xatom[c1], yatom[c1], zatom[c1]);  
         //      }

        fclose(InProcesFile);
        getatomcoor = 1;
        return getatomcoor;

}   /********************** end of getatomcoor ******************************/

/*****************************************************************************
this function gives the distances between two atoms in cartesian coordinates
*****************************************************************************/
double  atomsdistc(double xrotaxis1, double y1, double zrotaxis1, double x2, double y2, double z2)
{
    double distance;
        
    distance = sqrt((xrotaxis1-x2)*(xrotaxis1-x2) + (y1-y2)*(y1-y2) + (zrotaxis1-z2)*(zrotaxis1-z2));

    return distance;
}   /*************** end of atomscdist **********************/

/*****************************************************************************
             minor changes on the original function by Aileen Grey

this function gives the distances between two atoms of a periodic system
the function needs the lattice parameters given by an external variable
 extern double cell[N];

*****************************************************************************/
double  atomsdist(double x1, double y1, double z1, double x2, double y2, double z2)
{
        double atom[N];
        double o_atom[N];
        double per_atom[N];
        double ouratom[N];
        double o_ouratom[N];
        char line[201];
        char line2[201];
        char element[10];
        char type[10];
        double dist[10];
        double mag_dist;
        int i;
        int atom_num;

double  mag(double x, double y, double z);
double  coordtovec(int i, double cell[N], double per_atom[N]);
double  vectocoord(int i, double cell[N], double o_atom[N]);

 extern double cell[N];
        
        atom[1] = x1;
        atom[2] = y1;
        atom[3] = z1;
        ouratom[1] = x2;
        ouratom[2] = y2;
        ouratom[3] = z2;

        for (i=1;i<=3;i++)
          {
          if ((atom[i]-ouratom[i])>0.5)
             {
             per_atom[i]=(atom[i]-1.0);
             }
          else
             {
             if ((ouratom[i]-atom[i])>0.5)
               {
               per_atom[i]=(1.0+atom[i]);
               }
             else
               {
               per_atom[i]=atom[i];
               }
             }
          }
        for (i=1;i<=3;i++)
         {
         o_ouratom[i]=coordtovec(i, cell, ouratom);
         o_atom[i]=coordtovec(i, cell, per_atom);
         dist[i]=o_ouratom[i]-o_atom[i];
         }
         mag_dist = mag(dist[1], dist[2], dist[3]);

         return mag_dist;
}   /*************** end of the main body of atomsdist **********************/

/* function to find the magnitude of a vector */ 
double mag(double x, double y, double z)
{ 
          double mag;
          mag = sqrt((x*x)+(y*y)+(z*z)); 
          return mag;
}

/* function to calculate the cartesian coordinates relative to a */
/* orthonormal basis set with the x axis parallel to the x axis */
/* specified - units are in amstrongs */
double coordtovec(int i, double cell[N], double per_atom[N] )
{
  double a[N][N];
  double coordtovec;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));  
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3])); 

  coordtovec = (a[i][1]*per_atom[1])+(a[i][2]*per_atom[2])+(a[i][3]*per_atom[3]);

  return (coordtovec);
}

/*and function to change them back again */
double vectocoord(int i, double cell[N], double o_atom[N])
{
  double a[N][N];
  double acof[N][N];
  double ainv[N][N];
  double deta;
  double vectocoord;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3]));

  acof[1][1] = (a[2][2]*a[3][3])-(a[2][3]*a[3][2]);
  acof[1][2] = (a[1][3]*a[3][2])-(a[1][2]*a[3][3]);
  acof[1][3] = (a[1][2]*a[2][3])-(a[1][3]*a[2][2]);
  acof[2][1] = (a[2][3]*a[3][1])-(a[2][1]*a[3][3]);
  acof[2][2] = (a[1][1]*a[3][3])-(a[1][3]*a[3][1]);
  acof[2][3] = (a[1][3]*a[2][1])-(a[1][1]*a[2][3]);
  acof[3][1] = (a[2][1]*a[3][2])-(a[2][2]*a[3][1]);
  acof[3][2] = (a[1][2]*a[3][1])-(a[1][1]*a[3][2]);
  acof[3][3] = (a[1][1]*a[2][2])-(a[1][2]*a[2][1]);

  deta = (a[1][1]*acof[1][1])+(a[1][2]*acof[2][1])+(a[1][3]*acof[3][1]);

  ainv[1][1] = (1/deta)*(acof[1][1]);
  ainv[1][2] = (1/deta)*(acof[1][2]);
  ainv[1][3] = (1/deta)*(acof[1][3]);
  ainv[2][1] = (1/deta)*(acof[2][1]);
  ainv[2][2] = (1/deta)*(acof[2][2]);
  ainv[2][3] = (1/deta)*(acof[2][3]);
  ainv[3][1] = (1/deta)*(acof[3][1]);
  ainv[3][2] = (1/deta)*(acof[3][2]);
  ainv[3][3] = (1/deta)*(acof[3][3]);

  vectocoord = (ainv[i][1]*o_atom[1])+(ainv[i][2]*o_atom[2])+(ainv[i][3]*o_atom[3]);


  return (vectocoord);

}
/************************** end of atomsdist ********************************/

/*****************************************************************************
function gives 0 if the line not countain words, 1 else
*****************************************************************************/
int    blanckline(char line[Maxline])
{
        int    blanckline;
        int    c1, c2;
        
        blanckline = 0;
/*        printf(" in blancklin0 input:%s\n", line);     */
        if (strlen(line) == 1)
        {
/*           printf(" in blancklin1\n");                 */
           if (isspace(line[0]) == 0)
           {
             blanckline = 1;
           }
           else
           {
             blanckline = 0;
           }
        }
        else
        {
          for (c1 = 0; c1 <= strlen(line)-1;c1++)
            {
               if (isspace(line[c1]) == 0)
                 {
                   blanckline = 1;
                   c1 = strlen(line)+1;
                 }
            }       
        }


        return blanckline;

}    /***********************  end of blanckline  ****************************/

/*****************************************************************************
 function to calculate the cartesian coordinates relative to a 
 orthonormal basis set with the x axis parallel to the x axis 
 specified - units are in amstrongs 
 i = 1, 2 and 3 stand for  x, y and z, respectively
 
 original by Aileen Grey
*****************************************************************************/ 
double Cryst2Cartes(int i, double cell[N], double per_atom[N] )
{
  double a[N][N];
  double coordtovec;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));  
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3])); 

  coordtovec = (a[i][1]*per_atom[1])+(a[i][2]*per_atom[2])+(a[i][3]*per_atom[3]);

  return (coordtovec);
}  /**************  end of Cryst2Cartes ********************/ 

/*****************************************************************************
 function to calculate the crystallographic coordinates from cartesians
 i = 1, 2 and 3 stand for  x, y and z, respectively

 original by Aileen Grey
*****************************************************************************/
double Cartes2Cryst(int i, double cell[N], double per_atom[N] )
{
  double a[N][N];
  double acof[N][N];
  double ainv[N][N];
  double deta;
  double vectocoord;

  a[1][1]=cell[1];
  a[1][2]=cell[2]*cos(cell[6]*M_PI/180.0);
  a[1][3]=cell[3]*cos(cell[5]*M_PI/180.0);
  a[2][1]= 0.0;
  a[2][2]=cell[2]*sin(cell[6]*M_PI/180.0);
  a[2][3]=cell[3]*(cos(cell[4]*M_PI/180.0)
    -(cos(cell[5]*M_PI/180.0)*cos(cell[6]*M_PI/180.0)))
    /(sin(cell[6]*M_PI/180.0));
  a[3][1]=0.0;
  a[3][2]=0.0;
  a[3][3]=sqrt(cell[3]*cell[3]-(a[1][3]*a[1][3])-(a[2][3]*a[2][3]));

  acof[1][1] = (a[2][2]*a[3][3])-(a[2][3]*a[3][2]);
  acof[1][2] = (a[1][3]*a[3][2])-(a[1][2]*a[3][3]);
  acof[1][3] = (a[1][2]*a[2][3])-(a[1][3]*a[2][2]);
  acof[2][1] = (a[2][3]*a[3][1])-(a[2][1]*a[3][3]);
  acof[2][2] = (a[1][1]*a[3][3])-(a[1][3]*a[3][1]);
  acof[2][3] = (a[1][3]*a[2][1])-(a[1][1]*a[2][3]);
  acof[3][1] = (a[2][1]*a[3][2])-(a[2][2]*a[3][1]);
  acof[3][2] = (a[1][2]*a[3][1])-(a[1][1]*a[3][2]);
  acof[3][3] = (a[1][1]*a[2][2])-(a[1][2]*a[2][1]);

  deta = (a[1][1]*acof[1][1])+(a[1][2]*acof[2][1])+(a[1][3]*acof[3][1]);

  ainv[1][1] = (1/deta)*(acof[1][1]);
  ainv[1][2] = (1/deta)*(acof[1][2]);
  ainv[1][3] = (1/deta)*(acof[1][3]);
  ainv[2][1] = (1/deta)*(acof[2][1]);
  ainv[2][2] = (1/deta)*(acof[2][2]);
  ainv[2][3] = (1/deta)*(acof[2][3]);
  ainv[3][1] = (1/deta)*(acof[3][1]);
  ainv[3][2] = (1/deta)*(acof[3][2]);
  ainv[3][3] = (1/deta)*(acof[3][3]);

  vectocoord = (ainv[i][1]*per_atom[1])+(ainv[i][2]*per_atom[2])+(ainv[i][3]*per_atom[3]);


  return (vectocoord);

}  /***************** end of Cartes2Cryst ************************/    



/****************************************************************************
simple function to gives absolute value of x
****************************************************************************/
double   absx(double x)
{
  double x1;
  
  if (x >= 0)
  {
    x1 = x;
  }
  else
  {
    x1 = -x;
  }
 
  return x1;
  
}  /***************** end of absx ************************/ 

/*****************************************************************************
this function gives the angles between three atoms of a periodic system
the function needs the lattice parameters given by an external variable
extern double cell[N];

the formula is that by J. Buerger, appeared in his book "", chapter 23, pp. 629

x11, y11, z11 are the coordinates of the central atom of the angle
******************************************************************************/
double  atomsangle(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33)
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    double s12, s13;
    double deltax12, deltax13, deltay12, deltay13, deltaz12, deltaz13;
    double angle1;

    /**  first put all the atoms in the same unit cell, the next to the (0,0,0) one at the positive site **/
    if ( (x11 > 1.0) || (x11 < 0.0) )
    {
      if (x11 > 1.0)
      {
       x1 = x11 - 1.0;
      }
      else
      {
        x1 = x11 + 1.0;
      } 
    }
    else
    {
     x1 = x11;
    }
    
    if ( (x22 > 1.0) || (x22 < 0.0) )
    {
      if (x22 > 1.0)
      {
       x2 = x22 - 1.0;
      }
      else
      {
        x2 = x22 + 1.0;
      } 
    }
    else
    {
     x2 = x22;
    }
    
    if ( (x33 > 1.0) || (x33 < 0.0) )
    {
      if (x33 > 1.0)
      {
       x3 = x33 - 1.0;
      }
      else
      {
        x3 = x33 + 1.0;
      } 
    }
    else
    {
     x3 = x33;
    }
    
    if ( (y11 > 1.0) || (y11 < 0.0) )
    {
      if (y11 > 1.0)
      {
       y1 = y11 - 1.0;
      }
      else
      {
        y1 = y11 + 1.0;
      } 
    }
    else
    {
     y1 = y11;
    }
    
    if ( (y22 > 1.0) || (y22 < 0.0) )
    {
      if (y22 > 1.0)
      {
       y2 = y22 - 1.0;
      }
      else
      {
        y2 = y22 + 1.0;
      } 
    }
    else
    {
     y2 = y22;
    }
    
    if ( (y33 > 1.0) || (y33 < 0.0) )
    {
      if (y33 > 1.0)
      {
       y3 = y33 - 1.0;
      }
      else
      {
        y3 = y33 + 1.0;
      } 
    }
    else
    {
     y3 = y33;
    }
    
    if ( (z11 > 1.0) || (z11 < 0.0) )
    {
      if (z11 > 1.0)
      {
       z1 = z11 - 1.0;
      }
      else
      {
        z1 = z11 + 1.0;
      } 
    }
    else
    {
     z1 = z11;
    }
    
    if ( (z22 > 1.0) || (z22 < 0.0) )
    {
      if (z22 > 1.0)
      {
       z2 = z22 - 1.0;
      }
      else
      {
        z2 = z22 + 1.0;
      } 
    }
    else
    {
     z2 = z22;
    }
    
    if ( (z33 > 1.0) || (z33 < 0.0) )
    {
      if (z33 > 1.0)
      {
       z3 = z33 - 1.0;
      }
      else
      {
        z3 = z33 + 1.0;
      } 
    }
    else
    {
     z3 = z33;
    }
    /**  end of  put all the atoms in the same unit cell, the next to the (0,0,0) one at the positive site **/
    
    /**  first put the second atom coordinates as the closer distances of the first atom **/
    if ( (x1 - x2 > 0.5) || (x1 - x2 < -0.5) )
    {
      if (x1 - x2 > 0.5) 
      {
        x2 = x2 + 1.0;
      }
      else
      {
        x2 = x2 - 1.0;
      }
    }
    
    if ( (y1 - y2 > 0.5) || (y1 - y2 < -0.5) )
    {
      if (y1 - y2 > 0.5)
      {
        y2 = y2 + 1.0;
      }
      else
      {
        y2 = y2 - 1.0;
      }
    }
    
    if ( (z1 - z2 > 0.5) || (z1 - z2 < -0.5) )
    {
      if (z1 - z2 > 0.5)
      {
        z2 = z2 + 1.0;
      }
      else
      {
        z2 = z2 - 1.0;
      }
    }
    /**  end of put the second atom coordinates as the closer distances of the first atom **/
    
    /**  first put the third atom coordinates as the closer distances of the first atom **/
    if ( (x1 - x3 > 0.5) || (x1 - x3 < -0.5) )
    {
      if (x1 - x3 > 0.5) 
      {
        x3 = x3 + 1.0;
      }
      else
      {
        x3 = x3 - 1.0;
      }
    }
    
    if ( (y1 - y3 > 0.5) || (y1 - y3 < -0.5) )
    {
      if (y1 - y3 > 0.5)
      {
        y3 = y3 + 1.0;
      }
      else
      {
        y3 = y3 - 1.0;
      }
    }
    
    if ( (z1 - z3 > 0.5) || (z1 - z3 < -0.5) )
    {
      if (z1 - z3 > 0.5)
      {
        z3 = z3 + 1.0;
      }
      else
      {
        z3 = z3 - 1.0;
      }
    }

    /**  end of put the third atom coordinates as the closer distances of the first atom **/
    
    /**     this is the key segment of the function       **/
    s12 = atomsdist(x1, y1, z1, x2, y2, z2);
    s13 = atomsdist(x1, y1, z1, x3, y3, z3);    
    
/*    printf(" distance 1 %f    distance 2 %f\n", s12, s13); */
    
    deltax12 = x1 - x2;
    deltax13 = x1 - x3;
    deltay12 = y1 - y2;
    deltay13 = y1 - y3;
    deltaz12 = z1 - z2;
    deltaz13 = z1 - z3;
/*    
    printf(" deltax12 %f    deltax13  %f\n", deltax12, deltax13);
    printf(" deltay12 %f    deltay13  %f\n", deltay12, deltay13);
    printf(" deltaz12 %f    deltaz13  %f\n", deltaz12, deltaz13);
*/    
    angle1 = deltax12*deltax13*cell[1]*cell[1] + deltay12*deltay13*cell[2]*cell[2] + deltaz12*deltaz13*cell[3]*cell[3];
/*    printf(" 1            angle %f\n", angle1); */
    angle1 = angle1 + (deltax12*deltay13 + deltay12*deltax13)*cell[1]*cell[2]*cos(cell[6]*M_PI/180.0);
/*    printf(" 2            angle %f\n", angle1); */
    angle1 = angle1 + (deltaz12*deltax13 + deltax12*deltaz13)*cell[1]*cell[3]*cos(cell[5]*M_PI/180.0);
/*    printf(" 3            angle %f\n", angle1); */
    angle1 = angle1 + (deltay12*deltaz13 + deltaz12*deltay13)*cell[2]*cell[3]*cos(cell[4]*M_PI/180.0);
/*    printf(" 4            angle %f\n", angle1); */
    
    angle1 = angle1/(s12*s13);
/*    printf(" cosin of the angle %f\n", angle1); */
    
    // correcting 180 and 0 degrees for error induced by numerical truncation 
    if ( (angle1<-1.0) &&  (angle1>-1.01) )
    {
      angle1 = -1.0;
    }
    else
    if ( (angle1>1.0) &&  (angle1<1.01) )
    {
      angle1 = 1.0;
    }
	angle1 = acos(angle1)*180.0/M_PI;
/*    printf("              angle %f\n", angle1); */
    
    
    return angle1;
}
/************************** end of atomsangle ********************************/


/*****************************************************************************
this function gives the angles between three atoms in a cartesian non-periodic
system
x11, y11, z11 are the coordinates of the central atom of the angle
******************************************************************************/
double  atomsanglec(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33)
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
    double s12, s13, r12xr13;
    double deltax12, deltax13, deltay12, deltay13, deltaz12, deltaz13;
    double angle1;

    s12 = atomsdistc(x11, y11, z11, x22, y22, z22);
    s13 = atomsdistc(x11, y11, z11, x33, y33, z33);    
    
/*    printf(" distance 1 %f    distance 2 %f\n", s12, s13); */
    
    r12xr13 = (x22-x11)*(x33-x11) + (y22-y11)*(y33-y11) + (z22-z11)*(z33-z11);

    angle1 = r12xr13/(s12*s13);
/*    printf(" cosin of the angle %f\n", angle1); */
    
    // correcting 180 and 0 degrees for error induced by numerical truncation 
    if ( (angle1<-1.0) &&  (angle1>-1.01) )
    {
      angle1 = -1.0;
    }
    else
    if ( (angle1>1.0) &&  (angle1<1.01) )
    {
      angle1 = 1.0;
    }
    angle1 = acos(angle1)*180.0/M_PI;
/*    printf("              angle %f\n", angle1); */
    
    
    return angle1;
}
/************************** end of atomsanglec ********************************/

/****************************************************************************
  this function gives the nth word in a string
*****************************************************************************/
char  *grepnword(char *line1, int position_line, int nthword)
{
        int    c1, c2, c3, tempcount;
        char  grepnword1[Maxline];

        c2 = position_line;
        
        for (tempcount = 1; tempcount <= nthword; tempcount++)
        {
          c1 = c2;
          while ( (line1[c1] == ' ' )  && (c1 <= (strlen(line1)-1) ) )
             {
             c1++;
             }
          c2 = c1+1;
          while ( (line1[c2] != ' ')  && (c2 <= (strlen(line1)-1) ) )
             {
               c2++;
             }
          if (c2 == (strlen(line1)))
             {
               if (tempcount < nthword)
               {
                 tempcount = nthword+2;
                 strcpy(grepnword1, "");
               }
             }
        }

        if (tempcount != nthword+2)
        {
          for (c3 = 0; c3 < c1; c3++)
            line1++;
          sscanf(line1,"%s", grepnword1);
        }

        return grepnword1;

} /***********************  end of grepnword *******************************/

/*****************************************************************************
function to grep a string (line2) in a file,
returns 0  if the string is in the file, and
returns -1 if the string isn't in the file  
also returns as 'line' the line read, if exits !
*****************************************************************************/
int    strgrepf(char line1[400], char line2[400], int maxline1, FILE *filename1) 
{   
    char  *chtmp1;
    int   auxcount1;
    int   strgrepf;
    
    auxcount1 = 0;
    chtmp1 = fgets(line1, maxline1, filename1);
    while ((strstr(line1, line2) == NULL) && (chtmp1 != NULL))
        {                 
          chtmp1 = fgets(line1, maxline1, filename1);
        }
    if (chtmp1 != NULL)    
       strgrepf = 0;
    else 
       strgrepf = -1;

   return strgrepf;            
}    /***********************  end of strgrepf  ****************************/

/****************************************************************************
     this function gives the number of words in a line between
     two given points; count1pos, count2pos.
*****************************************************************************/
int    countwords(char line1[Maxline], int count1pos, int count2pos)
{
        int   c1, state_word;
        int   countwords;

        
        if ( isspace(line1[count1pos]) == 0) 
          {
            state_word = 1;
            countwords = 1;
          }
          else
          {
            state_word = 0;
            countwords = 0;
          }
          
          c1 = count1pos-1;
          while (c1 < count2pos)
            {
              c1++;
              if (isspace(line1[c1]) == 0)
                {
                  if (state_word == 0)
                    {
                      state_word = 1;
                      countwords++;
                    }
                }
                else
                {
                  if (state_word != 0)
                    state_word = 0;
                }   
            }

        return countwords;

}   /********************** end of countwords ******************************/

/****************************************************************************
  this function gives the word located between position n and m of a string
*****************************************************************************/
char  *grepnmword(char *line1, int position_ini, int position_final)
{
        int    c1, c2, c3, tempcont;
        char  grepnmword1[Maxline];
        char  grepnmword2[Maxline];
        char  grepnmword3[Maxline];
        char  ch1;
        int   ch_int1, ch_int2;

        strcpy(grepnmword2, "");

        if ((position_ini > strlen(line1)) || (position_ini > position_final) )
        {
          strcpy(grepnmword2, "");
        }
        else
        {
          if (position_final > strlen(line1))
            position_final = strlen(line1);
          for (c3 = 0; c3 < position_ini; c3++)
            line1++;

          ch_int1 = line1[0];
          sprintf(grepnmword2, "%c", ch_int1);
          sprintf(grepnmword3, "%c", ch_int1);



          c2 = position_final - position_ini;
          for (c1 = 1; c1 < c2; c1++)
          {
            line1++;
            ch_int1 = line1[0];
            sprintf(grepnmword3, "%c", ch_int1);
            strcat(grepnmword2, grepnmword3);
          }
        }

        return grepnmword2;

} /***********************  end of grepnmword *******************************/
                                                                               
