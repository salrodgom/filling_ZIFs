/****************************************************************************

-----------------------------------------------------------------------------
--------------------           zif_dist_angle_01         -------------------
----              Univ. Pablo Olavide, RASPA Group, Spain                ----
--------------------------        May 2013       ----------------------------                              
-----------------------------------------------------------------------------

  utilitary program to extract ZnN distances, NZnN angles, ZnNZn angles and Zn-Zn
  distances for a ZIF 
 
  crash_label is an error message flag: = 0 if error ocurrs, > 0 else 
*****************************************************************************/

# include <math.h>
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>

# define MaxAtomNum   10000
# define MaxAtomNum2   5000
# define MaxAtomNum3   5000
# define N               10
# define Maxline        200

FILE  *InpF;
FILE  *OutputFile;
char   OutputFName[Maxline], CoordFName[Maxline];
double xatom[MaxAtomNum], yatom[MaxAtomNum], zatom[MaxAtomNum], atomcharge[MaxAtomNum];
char   atomname[MaxAtomNum][8];
char   atom_c_s[MaxAtomNum];
int    atomnumber, CoorFType;
double cell[N];
double ZnNdistmin, ZnNdistmax;
double ZnNdistave, NZnNangleave, ZnNC1angleave, ZnNC2angleave;
int    countsingleT, countN;
int    singleTatoms[5][MaxAtomNum2];
int    N_neighb[6][MaxAtomNum3];
double imiddistmax, imiddistmin;
int    atomtype[MaxAtomNum], N_used[MaxAtomNum3];
// double NZnNanglemin, NZnNanglemax, ZnZndistmin, ZnZndistmax, ZnNZnanglemin, ZnNZnanglemax;
double deltaZnNdistave, deltaNZnNangleave, deltaZnNC1angleave, deltaZnNC2angleave;
double ZnNdistave1, NZnNangleave1, ZnNC1angleave1, ZnNC2angleave1;
double weighZnNdistave, weighNZnNangleave, weighZnNC1angleave, weighZnNC2angleave;

void   Initialization();
void   MainJob();
int    getatomcoor(char coordFName[Maxline], int CoorFType);
int    countwords(char line1[Maxline], int count1pos, int count2pos);
char  *grepnword(char *line1, int position_line, int nthword);
char  *grepnmword(char *line1, int position_ini, int position_final);
double atomsdist (double x1, double y1, double z1, double x2, double y2, double z2);
double atomsdistc(double x1, double y1, double z1, double x2, double y2, double z2);
double atomsangle(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33);
double atomsanglec(double x11, double y11, double z11, double x22, double y22, double z22, double x33, double y33, double z33);
double absx(double x);
void   FindSingleTzif();

void main(argc, argv)
        int argc;
        char **argv;
{
int    crash_label, count1, count2, count3;
double dist, angle_degree;
double cryst_atom[4];  

    crash_label = 0.0;
    
    InpF = fopen(argv[1],"r");
   
    Initialization();
    
    FindSingleTzif();
    
    MainJob(); 
    
    printf(" thanks for using zif_dist_angle_01 \n");
    
}     /******************************  end of main  **************************/

/******************************************************************************
function to initiate the job; reading general data and initialising
general variables.

******************************************************************************/
void    Initialization()
{   
    int  maxline, count1;
    char line[Maxline];

    maxline = Maxline;     
    
    fgets(line, maxline, InpF);
    sscanf(line,"%s %d", CoordFName, &CoorFType);
    /* name and type (gin, gout, xtl) of the input structure file */
    
    fgets(line, maxline, InpF);
    sscanf(line,"%s", OutputFName);
    /* name of Output file */
     
    fgets(line, maxline, InpF);
    sscanf(line,"%lf %lf", &ZnNdistmin, &ZnNdistmax); 
    /* min & max Zn-N distances */

    fgets(line, maxline, InpF);
    sscanf(line,"%lf %lf", &imiddistmin, &imiddistmax); 
    /* min & max imidazole intramolecular distances */     
    
    fclose(InpF);

    printf(" ms1 %s %d\n", CoordFName, CoorFType);
    printf(" ms2 %s %d\n", OutputFName);
    printf(" ms3 %f %f\n", ZnNdistmin, ZnNdistmax);
    printf(" ms4 %f %f\n", imiddistmin, imiddistmax);
    
    printf("\n");
    
}    /*  end of Initialization   */

/******************************************************************************
function to do the Main Job

******************************************************************************/
void    MainJob()
{
    int    maxline, Maincount1, count1, count2, count3, count4, count5, count6, temp_char;
    char   line[Maxline], auxs1[Maxline], auxs2[Maxline], auxs3[Maxline], ginfile[Maxline];
    int    num_ZnNdist, num_NZnNangle, num_ZnNC1angle, num_ZnNC2angle;
    double dist1, dist2, dist3, angle1, angle2;
    double value_runmin, value_runmax;
    double goodnessZnN, goodnessNZnN, goodnessZnNC1, goodnessZnNC2, goodnessoverall;
    
    maxline = Maxline;
    printf("Entering Distances and Angles Analysis\n");

    ZnNdistave1     = 2.0219; NZnNangleave1     = 109.429; ZnNC1angleave1     = 128.362; ZnNC2angleave1     = 125.867;
    deltaZnNdistave = 0.0316; deltaNZnNangleave = 5.30670; deltaZnNC1angleave = 9.39482; deltaZnNC2angleave = 9.28607;
    weighZnNdistave = 1.0   ; weighNZnNangleave = 0.3    ; weighZnNC1angleave = 0.7    ; weighZnNC2angleave = 0.7    ;   
 
    OutputFile = fopen(OutputFName,"w");
    
    // analyzing ZnN distances
    num_ZnNdist = 0; ZnNdistave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count1 = singleTatoms[0][Maincount1];
      printf("atom %5d ZnN distances:", count1 + 1);
      for (count2 = 1; count2 <= 4; count2++)
      {
        count3 = singleTatoms[count2][Maincount1];
        if (count3 > -1)
        {
          dist1 = atomsdist(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3]);
          num_ZnNdist++;
          ZnNdistave = ZnNdistave + dist1;
          printf(" %f", dist1);
          if (dist1<value_runmin)
          {
            value_runmin = dist1;
          }
          if (dist1>value_runmax)
          {
            value_runmax = dist1;
          }
        }
      }
      printf("\n");  
    }  
    ZnNdistave = ZnNdistave / num_ZnNdist;
    if ( (value_runmax - value_runmin) >= deltaZnNdistave )
    {
      goodnessZnN = (absx(ZnNdistave1 - ZnNdistave)/ZnNdistave1 + absx(deltaZnNdistave - value_runmax + value_runmin)/deltaZnNdistave) * weighZnNdistave;
    }
    else
    {
      goodnessZnN = (absx(ZnNdistave1 - ZnNdistave)/ZnNdistave1) * weighZnNdistave;
    }
    fprintf(OutputFile, "ZnN_distances %f %f %f goodnessZnN %f\n", value_runmin, ZnNdistave, value_runmax, goodnessZnN);
    // end of analyzing ZnN distances
    
    // analyzing NZnN angles
    num_NZnNangle = 0; NZnNangleave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countsingleT; Maincount1++)
    {
      count1 = singleTatoms[0][Maincount1];
      printf("atom %5d NZnN angles:", count1 + 1);
      for (count2 = 1; count2 < 4; count2++)
      {
        count3 = singleTatoms[count2][Maincount1];
        if (count3 > -1)
        {
          for (count4 = count2 + 1; count4 <= 4; count4++)
          { 
            count5 = singleTatoms[count4][Maincount1];
            if (count5 > -1)
            {
              angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
              num_NZnNangle++;
              NZnNangleave = NZnNangleave + angle1;
              printf(" %f", angle1);
              if (angle1<value_runmin)
              {
                value_runmin = angle1;
              }
              if (angle1>value_runmax)
              {
                value_runmax = angle1;
              }
            }
          }
        }
      }
      printf("\n");  
    }  
    NZnNangleave = NZnNangleave / num_NZnNangle;
    if ( (value_runmax - value_runmin) >= deltaNZnNangleave )
    {
      goodnessNZnN = (absx(NZnNangleave1 - NZnNangleave)/NZnNangleave1 + absx(deltaNZnNangleave - value_runmax + value_runmin)/deltaNZnNangleave) * weighNZnNangleave;
    }
    else
    {
      goodnessNZnN = (absx(NZnNangleave1 - NZnNangleave)/NZnNangleave1) * weighNZnNangleave;
    }
    fprintf(OutputFile, "NZnN_angles %f %f %f goodnessNZnN %f\n", value_runmin, NZnNangleave, value_runmax, goodnessNZnN);
    // end of analyzing NZnN angles
    
    // analyzing ZnNC1 angles
    num_ZnNC1angle = 0; ZnNC1angleave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countN; Maincount1++)
    {
      count1 = N_neighb[0][Maincount1];
      count3 = N_neighb[1][Maincount1];
      count5 = N_neighb[2][Maincount1];
      if ( (count5 > -1) && (count3 > -1) )
      {
        angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        num_ZnNC1angle++;
        ZnNC1angleave = ZnNC1angleave + angle1;
        printf("atom %5d ZnNC1 angle: %f\n", angle1);
        if (angle1<value_runmin)
        {
          value_runmin = angle1;
        }
        if (angle1>value_runmax)
        {
          value_runmax = angle1;
        }
      }
    }  
    ZnNC1angleave = ZnNC1angleave / num_ZnNC1angle;
    if ( (value_runmax - value_runmin) >= deltaZnNC1angleave )
    {
      goodnessZnNC1 = (absx(ZnNC1angleave1 - ZnNC1angleave)/ZnNC1angleave1 + absx(deltaZnNC1angleave - value_runmax + value_runmin)/deltaZnNC1angleave) * weighZnNC1angleave;
    }
    else
    {
      goodnessZnNC1 = (absx(ZnNC1angleave1 - ZnNC1angleave)/ZnNC1angleave1) * weighZnNC1angleave;
    }
    fprintf(OutputFile, "ZnNC1_angles %f %f %f goodnessZnNC1 %f \n", value_runmin, ZnNC1angleave, value_runmax, goodnessZnNC1);
    // end of analyzing ZnNC1 angles
    
    // analyzing ZnNC2 angles
    num_ZnNC2angle = 0; ZnNC2angleave = 0.0; value_runmin = 10000.0; value_runmax = -10000.0;
    for (Maincount1 = 0; Maincount1 <= countN; Maincount1++)
    {
      count1 = N_neighb[0][Maincount1];
      count3 = N_neighb[1][Maincount1];
      count5 = N_neighb[3][Maincount1];
      if ( (count5 > -1) && (count3 > -1) )
      {
        angle1 = atomsangle(xatom[count1], yatom[count1], zatom[count1], xatom[count3], yatom[count3], zatom[count3], xatom[count5], yatom[count5], zatom[count5]);
        num_ZnNC2angle++;
        ZnNC2angleave = ZnNC2angleave + angle1;
        printf("atom %5d ZnNC2 angle: %f\n", angle1);
        if (angle1<value_runmin)
        {
          value_runmin = angle1;
        }
        if (angle1>value_runmax)
        {
          value_runmax = angle1;
        }
      }
    }  
    ZnNC2angleave = ZnNC2angleave / num_ZnNC2angle;
    if ( (value_runmax - value_runmin) >= deltaZnNC2angleave )
    {
      goodnessZnNC2 = (absx(ZnNC2angleave1 - ZnNC2angleave)/ZnNC2angleave1 + absx(deltaZnNC2angleave - value_runmax + value_runmin)/deltaZnNC2angleave) * weighZnNC2angleave;
    }
    else
    if ( (value_runmax - value_runmin) >= deltaZnNdistave )
    {
      goodnessZnNC2 = (absx(ZnNC2angleave1 - ZnNC2angleave)/ZnNC2angleave1) * weighZnNC2angleave;
    }
    fprintf(OutputFile, "ZnNC2_angles %f %f %f goodnessZnNC2  %f\n", value_runmin, ZnNC2angleave, value_runmax, goodnessZnNC2);
    // end of analyzing ZnNC2 angles
    
    goodnessoverall = goodnessZnN + goodnessNZnN + goodnessZnNC1 + goodnessZnNC2;
    fprintf(OutputFile, "Overall_goodness %f\n", goodnessoverall); 
    fclose(OutputFile);

}  /************* end of MainJob *************/

/****************************************************************************
     this function extractes the atoms names, core-shel specification
     and coordinates

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
        int    count1, count2, count3, count4, count5;
        char  *chtmp1;

        printf("already in getatomcoor 0, reading file %s of type %d \n", coordFName, CoorFType);
        
        InProcesFile = fopen(coordFName,"r");
        getatomcoor = 1; 
        maxline    = Maxline;
        atomnumber = 0;
        readinline = 1;

        memset(line, '\0', sizeof(line));
        memset(s1, '\0', sizeof(s1)); memset(s2, '\0', sizeof(s2)); memset(s3, '\0', sizeof(s3)); memset(s4, '\0', sizeof(s4));
        memset(s5, '\0', sizeof(s5)); memset(s6, '\0', sizeof(s6)); memset(s7, '\0', sizeof(s7)); memset(s8, '\0', sizeof(s8));
        
        printf("already in getatomcoor 1, reading file %s of type %d \n", coordFName, CoorFType);
        
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
           sscanf(line,"%s", s1);
           while (strncmp(s1, "_cell_length_a", 14)!= 0)
           {
//             printf(" checking reading a cell parameter  %s", line);
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           memset(s1, '\0', sizeof(s1)); memset(s2, '\0', sizeof(s2));
           sscanf(line,"%s %s", s1, s2);
           printf(" checking 1 reading a cell parameter  %s %s", s2, line);
//           count3 = strlen(s2) - 1;
//           count3 = strlen(s2);
//           count2 = count3;
//           for (count1 = 0; count1 < count3; count1++)
//           {
//             if (s2[count1]=='(')
//             {
//               count2 = count1;
//               count1 = count3;
//             }
//           }
//           strcpy(s1, grepnmword(s2, 0, count2));

           count3 = strlen(s2);
           count2 = count3 - 1;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           memset(s1, '\0', sizeof(s1));
           strcpy(s1, grepnmword(s2, 0, count2));
           printf(" checking 2 reading a cell parameter  %s %s\n", s1, s2);
//           sscanf(s1, "%lf", &cell[1]);  OJO esto necesita una revision detallada porque falla al leer s1, recuerda OJO
           sscanf(s2, "%lf", &cell[1]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1)); memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
//           count3 = strlen(s2) - 1;
           count3 = strlen(s2);
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           memset(s1, '\0', sizeof(s1));
           strcpy(s1, grepnmword(s2, 0, count2));
           printf(" checking 2 reading b cell parameter  %s %s\n", s1, s2);
           sscanf(s1, "%lf", &cell[2]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1)); memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
//           count3 = strlen(s2) - 1;
           count3 = strlen(s2);
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           memset(s1, '\0', sizeof(s1));
           strcpy(s1, grepnmword(s2, 0, count2));
           printf(" checking 2 reading c cell parameter  %s %s\n", s1, s2);
           sscanf(s1, "%lf", &cell[3]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1)); memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile); 
           sscanf(line,"%s %s", s1, s2);
//           count3 = strlen(s2) - 1;
           count3 = strlen(s2);
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           memset(s1, '\0', sizeof(s1));
           strcpy(s1, grepnmword(s2, 0, count2));
           printf(" checking 2 reading alpha cell parameter  %s %s\n", s1, s2);
           sscanf(s1, "%lf", &cell[4]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1)); memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
//           count3 = strlen(s2) - 1;
           count3 = strlen(s2);
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           memset(s1, '\0', sizeof(s1));
           strcpy(s1, grepnmword(s2, 0, count2));
           printf(" checking 2 reading beta  cell parameter  %s %s\n", s1, s2);
           sscanf(s1, "%lf", &cell[5]);
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1)); memset(s2, '\0', sizeof(s2));
           fgets(line, maxline, InProcesFile);
           sscanf(line,"%s %s", s1, s2);
//           count3 = strlen(s2) - 1;
           count3 = strlen(s2);
           count2 = count3;
           for (count1 = 0; count1 < count3; count1++)
           {
             if (s2[count1]=='(')
             {
               count2 = count1;
               count1 = count3;
             }
           }
           printf(" checking 2 reading gamma cell parameter  %s %s\n", s1, s2);
           sscanf(s2, "%lf", &cell[6]);
           printf("%f %f %f %f %f %f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);

           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           //while (strncmp(s1, "_atom_site_fract_z", 18)!= 0)
           while (strncmp(s1, "_atom_site_", 11)!= 0)
           {
             fgets(line, maxline, InProcesFile);
             sscanf(line,"%s", s1);
           }
           count4 = 0; count5 = 0;
           if (strncmp(s1, "_atom_site_label", 16) == 0)
           {
             count4 = 1; count5 = 1;
           }
           else
           if (strncmp(s1, "_atom_site_type_symbol", 22) == 0)
           {
             count4 = 1; count5 = 2;
           }
           memset(line, '\0', sizeof(line));
           memset(s1, '\0', sizeof(s1));
           chtmp1 = fgets(line, maxline, InProcesFile);
           sscanf(line,"%s", s1);
           if (strncmp(s1, "_atom_site_type_symbol", 22) == 0)
           {
             count4++; // count5 = 2;
           }
           else
           if (strncmp(s1, "_atom_site_label", 22) == 0)
           {
             count4++; // count5 = 2;
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
               if (count4==1)
               {
                 sscanf(line,"%s %lf %lf %lf %lf", atomname[atomnumber], &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);  
               }
               else
               {
                 if (count5==1)
                 {
                   sscanf(line,"%s %s %lf %lf %lf %lf", atomname[atomnumber], s2, &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);
                 }
                 else
                 {
                   sscanf(line,"%s %s %lf %lf %lf %lf", s2, atomname[atomnumber], &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);
                 }
               } 
               // sscanf(line,"%s %s %lf %lf %lf %lf", atomname[atomnumber], s1, &xatom[atomnumber], &yatom[atomnumber], &zatom[atomnumber], &atomcharge[atomnumber]);
               printf("%d %s %.6f %.6f %.6f %.6f\n", atomnumber, atomname[atomnumber], xatom[atomnumber], yatom[atomnumber], zatom[atomnumber], atomcharge[atomnumber]);
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

           printf("already in getatomcoor 2\n");

           break;
        }
      
         printf("atomnumber %d\n", atomnumber); 
         printf("cell  %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]); 
         for (c1 = 0; c1 <= atomnumber-1; c1++)
               {
                 printf("%s %c %.5f %.5f %.5f\n", atomname[c1], atom_c_s[c1], xatom[c1], yatom[c1], zatom[c1]);  
               }

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

atom1 is the centre of the angle

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
    
    angle1 = acos(angle1)*180.0/M_PI;
/*    printf("              angle %f\n", angle1); */
    
    
    return angle1;
}
/************************** end of atomsangle ********************************/


/*****************************************************************************
this function gives the angles between three atoms in a cartesian non-periodic
system

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

        memset(grepnmword1, '\0', sizeof(grepnmword1));
        memset(grepnmword2, '\0', sizeof(grepnmword2));
        memset(grepnmword3, '\0', sizeof(grepnmword3));
        //strcpy(grepnmword2, "");

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
                                                                               
/****************************************************************************
     this function uses getatomcoor to extract the coordinates and to create
     a list of the singleT (stored in singleTatoms)
*****************************************************************************/
void   FindSingleTzif()
{
    int    count1, count2, count3, count4, count5, count6, count7;
    int    tmpAtoms[N];
    double disttmp, disttmp1;
    char   s1[Maxline], s2[Maxline];

    printf(" ms8 %s %d\n", CoordFName, CoorFType);
    
    countsingleT = 0; countN = 0;
//    count3 = (MaxAtomNum / 3) + 1;
    for (count1 = 0; count1 < MaxAtomNum2; count1++)
    {
      for (count2 = 0; count2 < 5; count2++)
      {
        singleTatoms[count2][count1] = -1;
      }
    }
    for (count1 = 0; count1 < MaxAtomNum3; count1++)
    {
      N_used[count1] = 0;
      for (count2 = 0; count2 < 6; count2++)
      {
        N_neighb[count2][count1] = -1;
      }
    }
    for (count1 = 0; count1 < MaxAtomNum; count1++)
    {
      atomtype[count1] = -1;
    }
    
    
    printf(" already in FindSingleTzif 1, ready to process file %s of type %d\n", CoordFName, CoorFType);
    
    count1 = getatomcoor(CoordFName, CoorFType);
    
/* 
    printf(" atom 0 is %s %4.5f %4.5f %4.5f \n", atomname[0], xatom[0], yatom[0], zatom[0]);
    printf(" atom 1 is %s %4.5f %4.5f %4.5f \n", atomname[1], xatom[1], yatom[1], zatom[1]);
    printf(" atom 2 is %s %4.5f %4.5f %4.5f \n", atomname[2], xatom[2], yatom[2], zatom[2]);
    printf(" atom 3 is %s %4.5f %4.5f %4.5f \n", atomname[3], xatom[3], yatom[3], zatom[3]);
    printf(" atom 4 is %s %4.5f %4.5f %4.5f \n", atomname[4], xatom[4], yatom[4], zatom[4]);
    printf(" atomnumber %d \n", atomnumber);
*/    
    printf(" already in FindSingleTzif 2, getatomcoor done!\n");
    printf(" cell %.5f %.5f %.5f %.5f %.5f %.5f\n", cell[1], cell[2], cell[3], cell[4], cell[5], cell[6]);

    for (count1 = 0; count1 < atomnumber; count1++)  // finding for T (Zn, Co, Cd) atoms
    {
      memset(s1, '\0', sizeof(s1));
      strcpy(s1, atomname[count1]); 
      if ( (strncmp(s1, "Zn", 2) == 0) || (strncmp(s1, "zn", 2) == 0) || (strncmp(s1, "Co", 2) == 0) || (strncmp(s1, "co", 2) == 0) || (strncmp(s1, "Cd", 2) == 0) || (strncmp(s1, "cd", 2) == 0) )
      {
        tmpAtoms[0] = count1; 
        tmpAtoms[1] = -1;    tmpAtoms[2] = -1;
        tmpAtoms[3] = -1;    tmpAtoms[4] = -1;

        //printf("scaning Zn atom %d coord %.5f %.5f %.5f\n", count1, xatom[count1], yatom[count1], zatom[count1]);
        count3 = 0; 
        for (count2 = 0; count2 < atomnumber; count2++)
        {
          memset(s2, '\0', sizeof(s2));
          strcpy(s2, atomname[count2]);
          if ( (strncmp(s2, "N", 1) == 0) || (strncmp(s2, "n", 1) == 0) )
          { 
            disttmp = atomsdist(xatom[count1], yatom[count1], zatom[count1], xatom[count2], yatom[count2], zatom[count2]);
            if ((disttmp<=ZnNdistmax) && (disttmp>=ZnNdistmin)) 
            {
              count3++;
              tmpAtoms[count3] = count2;
            }
            //printf("scaning N  atom %d dist %.5f coord %.5f %.5f %.5f\n", count2, disttmp, xatom[count2], yatom[count2], zatom[count2]);
          }
        }
        printf(" current atom %d with %d bonds\n", count1, count3);
        if (count3 == 4)
        {
          countsingleT++;
          singleTatoms[0][countsingleT] = tmpAtoms[0];
          singleTatoms[1][countsingleT] = tmpAtoms[1];
          singleTatoms[2][countsingleT] = tmpAtoms[2];
          singleTatoms[3][countsingleT] = tmpAtoms[3];
          singleTatoms[4][countsingleT] = tmpAtoms[4];
          
          count2 = singleTatoms[0][countsingleT];
          atomtype[count2] = 1;
          count2 = singleTatoms[1][countsingleT];
          atomtype[count2] = 2;
          count2 = singleTatoms[2][countsingleT];
          atomtype[count2] = 2;
          count2 = singleTatoms[3][countsingleT];
          atomtype[count2] = 2;
          count2 = singleTatoms[4][countsingleT];
          atomtype[count2] = 2;
        
          printf(" atom 0 of singleTatoms %d is %d atomtype %d\n", countsingleT, singleTatoms[0][countsingleT], atomtype[singleTatoms[0][countsingleT]]);
          printf(" atom 1 of singleTatoms %d is %d atomtype %d\n", countsingleT, singleTatoms[1][countsingleT], atomtype[singleTatoms[1][countsingleT]]);
          printf(" atom 2 of singleTatoms %d is %d atomtype %d\n", countsingleT, singleTatoms[2][countsingleT], atomtype[singleTatoms[2][countsingleT]]);
          printf(" atom 3 of singleTatoms %d is %d atomtype %d\n", countsingleT, singleTatoms[3][countsingleT], atomtype[singleTatoms[3][countsingleT]]);
          printf(" atom 4 of singleTatoms %d is %d atomtype %d\n", countsingleT, singleTatoms[4][countsingleT], atomtype[singleTatoms[4][countsingleT]]);
        }
      }
    }          // end of finding for T (Zn, Co, Cd) atoms

    for (count1 = 1; count1 <= countsingleT; count1++)  // finding for N atoms
    {
      // 0 = N itself, 1 = Zn, 2 = C1, 3 = C2, 4 = N other next to 2, 5 = C2 next to 4
      for (count2 = 1; count2 <= 4; count2++)  
      {
        countN++;

        N_neighb[0][countN] = singleTatoms[count2][count1];
        N_neighb[1][countN] = singleTatoms[0][count1];
        count5              = N_neighb[0][countN];
        
        tmpAtoms[0] = N_neighb[0][countN];  tmpAtoms[1] = -1;  tmpAtoms[2] = -1;
        count4 = 0;
        for (count3 = 0; count3 < atomnumber; count3++)  
        {
          strcpy(s1, atomname[count3]);
          if ( ( (strncmp(s1, "C", 1) == 0) || (strncmp(s1, "c", 1) == 0) ) && ( (strncmp(s1, "Co", 2) != 0) && (strncmp(s1, "co", 2) != 0) && (strncmp(s1, "Cd", 2) != 0) && (strncmp(s1, "cd", 2) != 0) ) )
          {
            disttmp = atomsdist(xatom[count5], yatom[count5], zatom[count5], xatom[count3], yatom[count3], zatom[count3]);
        //    printf(" checking 1 N_neighb Zn %d N %d C running %d name %s dist %.5f distcuoff %.5f %.5f\n", count1, count2, count3, s1, disttmp, imiddistmin, imiddistmax);
            if ((disttmp<=imiddistmax) && (disttmp>=imiddistmin)) 
            {
        //      printf(" checking 2 N_neighb Zn %d N %d C running %d name %s dist %.5f\n", count1, count2, count3, s1, disttmp);
              count4++;
              tmpAtoms[count4] = count3;
              if (count4 == 2)
              {
                count3 = atomnumber + 1;
              }
            }
          }
        }        

        count4 = tmpAtoms[1]; 
        for (count3 = 0; count3 < atomnumber; count3++)  
        {
          strcpy(s1, atomname[count3]);
          if ( ( (strncmp(s1, "C", 1) == 0) || (strncmp(s1, "c", 1) == 0) || (strncmp(s1, "N", 1) == 0) || (strncmp(s1, "n", 1) == 0) ) && ( (strncmp(s1, "Co", 2) != 0) && (strncmp(s1, "co", 2) != 0) && (strncmp(s1, "Cd", 2) != 0) && (strncmp(s1, "cd", 2) != 0) ) )
//          if ( (strncmp(s1, "N", 1) == 0) || (strncmp(s1, "n", 1) == 0) )
		  {
            disttmp1 = atomsdist(xatom[count4], yatom[count4], zatom[count4], xatom[count3], yatom[count3], zatom[count3]);
            if ((disttmp1<=imiddistmax) && (disttmp1>=imiddistmin)) 
            {
              if ( (strncmp(s1, "N", 1) == 0) || (strncmp(s1, "n", 1) == 0) )
              {
//                for (count5 = 0; count5 < atomnumber; count5++)  
//                {    
//                  strcpy(s2, atomname[count1]);
//                  if ( (strncmp(s2, "Zn", 2) == 0) || (strncmp(s2, "zn", 2) == 0) || (strncmp(s1, "Co", 2) == 0) || (strncmp(s1, "co", 2) == 0) || (strncmp(s1, "Cd", 2) == 0) || (strncmp(s1, "cd", 2) == 0) )
//                  {
//                  	disttmp1 = atomsdist(xatom[count4], yatom[count4], zatom[count4], xatom[count3], yatom[count3], zatom[count3]);
//                   if ((disttmp1<=imiddistmax) && (disttmp1>=imiddistmin)) 
//                    {
//                      N_neighb[2][countN] = tmpAtoms[1];
//                      N_neighb[3][countN] = tmpAtoms[2];
//                      N_neighb[4][countN] = count3;
//                      count5 = atomnumber + 1;
//                      atomtype[tmpAtoms[1]] = 3;
//                      atomtype[tmpAtoms[2]] = 4;
//                    }
//                  }
//                }
                if (atomtype[count3] == 2)
                {
                  N_neighb[2][countN] = tmpAtoms[1];
                  N_neighb[3][countN] = tmpAtoms[2];
                  N_neighb[4][countN] = count3;
                  atomtype[tmpAtoms[1]] = 3;
                  atomtype[tmpAtoms[2]] = 4;
                }
                else
                {
                  N_neighb[3][countN] = tmpAtoms[1];
                  N_neighb[2][countN] = tmpAtoms[2];
                  atomtype[tmpAtoms[1]] = 4;
                  atomtype[tmpAtoms[2]] = 3;
                }
              }
              else
              {
              	N_neighb[3][countN] = tmpAtoms[1];
                N_neighb[2][countN] = tmpAtoms[2];
                atomtype[tmpAtoms[1]] = 4;
                atomtype[tmpAtoms[2]] = 3;
              }
        //      printf(" checking 3 N_neighb Zn %d N %d C_o_N running %d name %s dist %.5f\n", count1, count2, count3, s1, disttmp1);
              // now finding N_neighb[4][countN] if still needed
              if (N_neighb[4][countN] == -1)
              {
			    count6 = N_neighb[2][countN];
                for (count5 = 0; count5 < atomnumber; count5++)  
                {
                  if ( (atomtype[count5] == 2) && (count5 != N_neighb[0][countN]) )
                  {
                  	disttmp1 = atomsdist(xatom[count6], yatom[count6], zatom[count6], xatom[count5], yatom[count5], zatom[count5]);
                    if ((disttmp1<=imiddistmax) && (disttmp1>=imiddistmin))
                    {
                      N_neighb[4][countN] = count5;
                      count5 = atomnumber + 1;
                    }
                  }
                }
              }
              else
              {
                disttmp1 = -1000.0;
              }
              count6 = N_neighb[4][countN];
        //      printf(" checking 4 N_neighb Zn %d N %d N_2 running %d name %s dist %.5f\n", count1, count2, count6, atomname[count6], disttmp1);
              for (count5 = 0; count5 < atomnumber; count5++)  
              {
                strcpy(s2, atomname[count5]);
                if ( ( (strncmp(s2, "C", 1) == 0) || (strncmp(s2, "c", 1) == 0) ) && (atomtype[count5] != 3) )
                {
                  disttmp1 = atomsdist(xatom[count6], yatom[count6], zatom[count6], xatom[count5], yatom[count5], zatom[count5]);
        //          printf(" checking 5 N_neighb Zn %d N %d C2 running %d name %s dist %.5f\n", count1, count2, count3, s2, disttmp1);
                  if ((disttmp1<=imiddistmax) && (disttmp1>=imiddistmin))
                  {
                    N_neighb[5][countN] = count5;
                    atomtype[count5] = 4;
        //            printf(" checking 6 N_neighb Zn %d N %d C2 running %d name %s dist %.5f N_neighb[5] %d atomtype[c2] %d\n", count1, count2, count3, s2, disttmp1, N_neighb[5][countN], atomtype[count5]);
                    count5 = atomnumber + 1;
                  }	
                }
              }  
            }
          }
        } 
		printf(" N neighbor %d of Zn %d\n", count2, count1);
        printf(" atom 0 of N_neighb %d is %3d atomtype %d name %s \n", countN, N_neighb[0][countN], atomtype[N_neighb[0][countN]], atomname[N_neighb[0][countN]]);
        printf(" atom 1 of N_neighb %d is %3d atomtype %d name %s \n", countN, N_neighb[1][countN], atomtype[N_neighb[1][countN]], atomname[N_neighb[1][countN]]);
        printf(" atom 2 of N_neighb %d is %3d atomtype %d name %s \n", countN, N_neighb[2][countN], atomtype[N_neighb[2][countN]], atomname[N_neighb[2][countN]]);
        printf(" atom 3 of N_neighb %d is %3d atomtype %d name %s \n", countN, N_neighb[3][countN], atomtype[N_neighb[3][countN]], atomname[N_neighb[3][countN]]);
        printf(" atom 4 of N_neighb %d is %3d atomtype %d name %s \n", countN, N_neighb[4][countN], atomtype[N_neighb[4][countN]], atomname[N_neighb[4][countN]]);
        printf(" atom 5 of N_neighb %d is %3d atomtype %d name %s \n", countN, N_neighb[5][countN], atomtype[N_neighb[5][countN]], atomname[N_neighb[5][countN]]);    
      }
    }          // end of finding for N atoms    
 
}   /********************** end of FindSingleTzif ******************************/
