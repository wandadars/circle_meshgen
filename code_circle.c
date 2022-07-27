// generating square quadrilateral everywhere
// generate linear mesh and rotate it to create 2D mesh
// removing contraint on farfield mesh being square, no stretching1
// removing specification of delr, rather specifying N
// making new mesh semicirclular block: dec 12th 2k6. 
// making new mesh single circular block
// changing circle calculation ...making it angel based 
// impreviing value of PI from 3.1416 to 3.1415926535897932384626433832792 
// improving upon precision while writing .. making %lf to %1.16e
// again modyfing topology
// implemented stretching properly, matching spacing between blocks
// modified output_format() 
// modifying earlier mesh for cylinder shock diffraction calculations
/**********************************************************************/
/*  Code to generate TFI grid in square domain with circular hole */
/*  Author  : Manoj Kumar Parmar                      */
/*  Date    : Dec 12th, 2006                          */
/**********************************************************************/
    
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define PI 3.1415926535897932384626433832792
#define NO_OF_COLUMNS 5
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define NEXT(x)   (x==3 ? 0 : x+1)
#define PREV(x)   (x==0 ? 3 : x-1)
#define DOTPRODUCT(A,B) ((A.x)*(B.x)+(A.y)*(B.y)+(A.z)*(B.z))

typedef struct
{
    double  x, y, z ;
    long int vertexno ;
} vertex ;

typedef struct
{
    long int Imax ;
    vertex  *edgegrid ;
    int gridded ;   /** ==1 when TFI is done **/
    double  zi_stretching ;
    int iscircle ; /** ==1 when edge is circular section **/
    vertex  centre ;
    double  radius ;

} edge ;

typedef struct
{
    int Block_id ;
    long int Imax, Jmax ;

    vertex  corner[4] ;
    edge    boundary[4] ;
    vertex  **grid ;

    int gridded ; /** ==1 when TFI is done **/
    double  zi_stretching, eta_stretching ;

        /* MULTI-BLOCK FEATURES */
    double  angle1, angle2 ;
    int direction1, direction2 ;
    int shared_side[4] ;
    int neighbor[4] ;

} block ;

double A11, A12, A13, A22, A23, A33 ; /** contravariant terms **/
double a11, a12, a13, a22, a23, a33 ; /** covariant terms **/
vertex    a1, a2, a3 ; /** covariant vectors **/

/* GLOBAL VARIABLE DECLARATIONS */
block     *Blocks ;
int    Block_num=1 ; /* Not dealing with multi-blocks yet */
long int   gpts=10, gpts1=10 ;
double    *LinearMesh ;

/* Global variable for cylinder boundary */
long int    ncylpts  ;
vertex      *cylinder ;

/* PARAMETERS FOR DOMAIN */
double   D=0.01225, R=0.006125 ; /* Diameter of cicular hole */
double   D1=1.225, R1=0.6125 ; /* Diameter of outer cicular boundary */
double   D1_D=100.0 ; /* Ratio of inner and outer boundary */
double   angle=2*PI, angle1=0.0 ;
double   stretching=1.0 ;
long int Imax=10, Jmax=10 ; 
int  mode=0 ; /* 0 automatic,  !=0 manual */

/* DECLARATIONS FOR FUNCTIONS */
void    compute_Jmax(void) ;
double  find_stretching(double, double, long int, double) ;
double  find_stretching1(double, double, long int, double) ;
void    impose_symmetry(void) ;
void    initialize_flags(void) ;
void    initialize_meshsize(void) ;
void    generate_2Dmesh(void) ;
void    generate_linearmesh(void) ;
int generate_grid(int) ;
int generate_TFI_block(int) ;
int generate_TFI_edge(int, int) ;
void    output_format(void) ;
void    output_plot3D(void) ;
void    output_tecplot(void) ;
int read_inputs(int, char**) ;
void    testcase(void) ;

/* NOT IN USE ANYMORE */
int     laplace_updatexyz(int, int, int) ;
int     runlaplacian(int) ;
int     runlaplacian1(int) ;


/* MAIN PROGRAM BIGINS HERE */
int main(int argc, char *argv[])
{
    int     j,error ;
    char    char1 ;

        /* Read in command line inputs */
        printf("\n") ;
        error = read_inputs(argc,argv) ;
        if(error != 0)
    {
        printf("read_inputs() crashed \n") ;
        return(error) ;
    }

        /* malloc for default mesh in single block */
        printf("\n") ;
        printf("Initializing default mesh in single block \n") ;
    Block_num = 1 ;
        Blocks = (block *) malloc (Block_num*sizeof(block)) ;
    initialize_flags() ;
    initialize_meshsize() ;
        generate_grid(0) ;

        /* generate 1D mesh */
        printf("\n") ;
        printf("generating 1D mesh \n") ;
    generate_linearmesh() ;

        /* generate 2D mesh. Works only for single block */
        printf("\n") ;
        printf("generating 2D mesh \n") ;
        generate_2Dmesh() ;

        /* impose symmetry in 2D mesh. Works only for single block */
        printf("\n") ;
        printf("Imposing symmetry in 2D mesh \n") ;
    impose_symmetry() ;

        printf("\n") ;
        printf("Writing plot3D file \n") ;
    output_plot3D() ;
        printf("Writing tecplot file \n") ;
    output_tecplot() ;
        printf("Writing mesh file for input to conflu \n") ;
    output_format() ;

    return (0) ;
}

void compute_Jmax()
{
    int i, count=0, do_loop=1 ;
    double dtheta, delr, distance ;
        FILE *fp ;

        fp = fopen("delr.dat","w") ;

        dtheta = 2*PI/(Imax-1.0) ;

        delr = 0 ;
        distance = 0 ;
        count = 1 ;

        fprintf(fp,"%d %lf\n",count,delr) ;

        do_loop = 1 ;
        while( do_loop == 1 )
        {
          if(R1-R-distance > delr)
          {
            delr = (R+distance)*dtheta ; 
            if(distance+delr > R1-R ) 
            {
              delr = (R1-R)-distance ; 
              distance = R1-R ;
              do_loop = 0 ;
            }
            else
            {
              distance = distance + delr ;
            }
            count++ ;
            fprintf(fp,"%d %lf\n",count,delr) ;
          }
          else
          {
            distance = R1-R ;
            do_loop = 0 ;
          }
        }

        printf("\n") ;
        printf("Jmax found to be %d \n",count) ;
        Jmax = count ;

        fclose(fp) ;

    return ;
}

void generate_linearmesh()
{
    int i, count=0, do_loop=1, need_smoothening=0, stencil_size=10 ;
    double dtheta, delr, distance, zi, r1, r2 ;

    LinearMesh = (double *) malloc (Jmax*sizeof(double)) ;

        dtheta = 2*PI/(Imax-1.0) ;

        delr = 0 ;
        distance = 0 ;
        count = 1 ;
        LinearMesh[count-1] = R+distance ;

        do_loop = 1 ; 
        need_smoothening = 0 ;
        while( do_loop == 1 )
        {
          if(R1-R-distance > delr)
          {
            delr = (R+distance)*dtheta ; 
            if(distance+delr > R1-R ) 
            {
              distance = R1-R ;
              do_loop = 0 ;
            }
            else
            {
              distance = distance + delr ;
            }
            count++ ;
            LinearMesh[count-1] = R+distance ;
          }
          else
          {
            distance = R1-R ;
            LinearMesh[count-1] = R+distance ;
            need_smoothening = 1 ;
            do_loop = 0 ;
          }
        }

        /** smoothen the mesh at farfield **/
        if(need_smoothening == 1)
        {
           stencil_size = Jmax/10 ;
           if(stencil_size > 10 ) stencil_size = 10 ; 
           if(stencil_size < 3 ) stencil_size = 3 ; 
           if(stencil_size >= Jmax ) need_smoothening = 0 ; 
        }
        if(need_smoothening == 1)
        {
           r1 = LinearMesh[Jmax-stencil_size] ;
           r2 = LinearMesh[Jmax-1] ;

           for(i=Jmax-stencil_size; i<Jmax; i++)
           {
             zi = (i*1.0 - Jmax+stencil_size)/(stencil_size-1.0) ;
             LinearMesh[i] = r1*(1-zi) + r2*(zi) ;
           }
        }

    return ;
}

void generate_2Dmesh()
{
    int i,j ;
    double zi,theta ;

        /*** generate 2D domain by rotating LinearMesh **********/
    for(i=0; i<Imax; i++)
    {
        zi = (i*1.0)/(Imax-1.0) ;
        theta = (1.0-zi)*angle + (zi)*angle1 ;

        for(j=0; j<Jmax; j++)
        {
            Blocks[0].grid[i][j].x =  cos(theta)*LinearMesh[j] ;
            Blocks[0].grid[i][j].y =  sin(theta)*LinearMesh[j] ;
            Blocks[0].grid[i][j].z = 0.0 ;
        }
    }
    return ;
}

void impose_symmetry()
{
    int i,j ;

        /*** imposing y=0 on center line ********/
    if((Imax-1)%2 == 0)
    {
      for(j=0; j<Jmax; j++)
      {
        Blocks[0].grid[(Imax/2)][j].y = 0.0 ;
      }
    }

        /*** imposing y=0 on (i=0, i=Imax-1) ********/
        i = 0 ; 
    for(j=0; j<Jmax; j++) Blocks[0].grid[i][j].y = 0.0 ;
        i = Imax-1 ;    
    for(j=0; j<Jmax; j++) Blocks[0].grid[i][j].y = 0.0 ;

        /*** imposing symmetry in y>0 and y<0 regions  ********/
    for(i=0; i<Imax; i++)
    {
      for(j=0; j<Jmax; j++)
      {
        Blocks[0].grid[Imax-1-i][j].x =  Blocks[0].grid[i][j].x ;
        Blocks[0].grid[Imax-1-i][j].y = -Blocks[0].grid[i][j].y ;
        Blocks[0].grid[Imax-1-i][j].z =  Blocks[0].grid[i][j].z ;
      }
    }

        /*** imposing x=0 on quarter and 3 quarter line ********/
    if((Imax-1)%4 == 0)
    {
      for(j=0; j<Jmax; j++)
      {
        Blocks[0].grid[  (Imax/4)][j].x = 0.0 ;
        Blocks[0].grid[3*(Imax/4)][j].x = 0.0 ;
      }
    }

        /*** imposing symmetry in x>0 and x<0 regions ********/
    if((Imax-1)%4 == 0)
    {
      for(i=1; i<Imax/4; i++)
          {
        for(j=0; j<Jmax; j++)
        {
        Blocks[0].grid[  (Imax/4)-i][j].x = -Blocks[0].grid[  (Imax/4)+i][j].x ;
        Blocks[0].grid[  (Imax/4)-i][j].y =  Blocks[0].grid[  (Imax/4)+i][j].y ;
        Blocks[0].grid[3*(Imax/4)-i][j].x = -Blocks[0].grid[3*(Imax/4)+i][j].x ;
        Blocks[0].grid[3*(Imax/4)-i][j].y =  Blocks[0].grid[3*(Imax/4)+i][j].y ;
        }
          }
    }

    return ;
}

void initialize_meshsize()
{
    int i, j ;

    /** initializing No of grid points **/
        for(j=0; j<Block_num; j++)
    { 
        Blocks[j].Imax = Imax ; 
        Blocks[j].Jmax = Jmax ;
    }
        for(j=0; j<Block_num; j++)
    { 
        Blocks[j].boundary[0].Imax = Blocks[j].Imax ;
        Blocks[j].boundary[2].Imax = Blocks[j].Imax ;
            Blocks[j].boundary[1].Imax = Blocks[j].Jmax ;
        Blocks[j].boundary[3].Imax = Blocks[j].Jmax ;
    }

    return ;    
}

void initialize_flags()
{
    int i,j ;
    for(j=0; j<Block_num; j++)
    {
        Blocks[j].gridded = 0 ;
            for(i=0; i<4; i++) Blocks[j].boundary[i].gridded = 0 ;
    }
    return ;
}

/* Order in which input are read : mode, D, D1_D, Imax */
int read_inputs(int argc, char *argv[])
{

        switch(argc)
        {
        case 1 :
                printf("well i am taking default values for all parameters\n") ;
                printf("mode=%d,  D=%1.4e,  D1_D=%1.4e,  Imax=%d  \n",
                        mode,D,D1_D,gpts) ;
                break ;
        case 2 :
                mode = atoi(argv[1]) ;
                if(mode == 0)
                {
                 printf("you have provided just mode\n") ;
                 printf("mode=%d,  D=%1.4e,  D1_D=%1.4e,  Imax=%d  \n",
                        mode,D,D1_D,gpts) ;
                }
                break ;
        case 3 :
                mode = atoi(argv[1]) ;
                D = atof(argv[2]) ;
                if(mode == 0)
                {
                  printf("you have provided just mode and D\n") ;
                  printf("mode=%d,  D=%1.4e,  D1_D=%1.4e,  Imax=%d  \n",
                        mode,D,D1_D,gpts) ;
                }
                break ;
        case 4 :
                mode = atoi(argv[1]) ;
                D = atof(argv[2]) ;
                D1_D = atof(argv[3]) ;
                if(mode == 0)
                {
                  printf("you have provided mode, D and D1_D\n") ;
                  printf("mode=%d,  D=%1.4e,  D1_D=%1.4e,  Imax=%d  \n",
                        mode,D,D1_D,gpts) ;
                }
                break ;
        case 5 :
                mode = atoi(argv[1]) ;
                D = atof(argv[2]) ;
                D1_D = atof(argv[3]) ;
                gpts = atoi(argv[4]) ;
                if(mode == 0)
                {
                  printf("you have provided mode, D, D1_D and Imax\n") ;
                  printf("mode=%d,  D=%1.4e,  D1_D=%1.4e,  Imax=%d  \n",
                        mode,D,D1_D,gpts) ;
                }
                break ;
        default :
                printf("you are providing mode that 4 inputs\n") ;
                printf("you should provide \"mode\" \"D\" \"D1_D\" and \"Imax\"\n") ;
                return(1) ;
                break ;
        }

        /*** setting up domain **/
        R = D/2 ;
        D1 = D1_D*D ; R1 = D1/2 ;
        /************************/

        if(mode != 0)
        {
                printf("Mesh would be generated according parameters given by user\n") ;
                printf("All quadrilaterals may not be square depending on user given parameter \n") ;
                printf("enter \"Imax\" \"Jmax\" and \"stretching\" : ");
                scanf("%d %d %lf",&gpts,&gpts1,&stretching) ;
        }

        if(gpts < 2 || gpts1 < 2)
        {
                printf("need to increase resolution Imax=%d or Kmax=%d too low\n",gpts,gpts1) ;
                return(2) ;
        }
        Imax = gpts ;
        Jmax = gpts1 ;

        if(mode == 0)
        {
          printf("Jmax would be computed such that all quadrilaterals are square \n") ;
          printf("between inner radius %lf and outer radius %lf \n",R,R1) ;
          compute_Jmax() ;
        }

        printf("\n") ;
        printf("Imax=%d Jmax=%d \n",Imax,Jmax) ;
        if(mode != 0) printf("stretching=%lf \n",stretching) ;

        return(0) ;
}

double  find_stretching(double x1, double x2 , long int num, double delx)
{
    double ratio = 2.0 ;

    ratio = (log((fabs(x2-x1))/delx))/(log(num-1.0)) ;

//  if(1==2)
    {
    printf("x1=%lf x2=%lf num=%d delx=%lf\n",x1,x2,num,delx) ;
    printf("log(fabs(x2-x1)/delx)=%lf \n",log((fabs(x2-x1))/fabs(delx))) ;
    printf("log(num-1)=%lf \n",log(num-1)) ;
    printf("num=%d ratio = %lf\n",num,ratio) ;
    }

    return(ratio) ;
}

double  find_stretching1(double x1, double x2 , long int num, double delx)
{
        double ratio = 2.0 ;

        ratio = (log(1-(delx/fabs(x2-x1))))/(log((num-2.0)/(num-1.0))) ;

//      if(1==2)
        {
        printf("x1=%lf x2=%lf num=%d delx=%lf\n",x1,x2,num,delx) ;
        printf("log(1-delx/fabs(x2-x1))=%lf \n",log(1-fabs(delx)/fabs(x2-x1))) ;
        printf("log((num-2)/(num-1))=%lf \n",log((num-2.0)/(num-1.0))) ;
        printf("num=%d ratio = %lf\n",num,ratio) ;
        }

        return(ratio) ;
}

int generate_grid(int blockno)
{
    int i, error=0 ;
 
        for(i=0; i<4; i++)
    { 
        error = generate_TFI_edge(blockno,i) ;
        if(error != 0)
        {
            printf("generate_TFI_edge() crashed\n") ;
            return(error) ;
        }
    }
    error = generate_TFI_block(blockno) ;
    if(error != 0)
    {
        printf("generate_TFI_block() crashed\n") ;
        return(error) ;
    }

    return(error) ;
}       

int generate_TFI_edge(int block, int side)
{
    long int i, M ;
    double  zi, x1, x2, y1, y2, z1, z2 ;
    double  zi_stretching = 1 ;

    if(Blocks[block].boundary[side].gridded == 1) return(0) ;

//  free(Blocks[block].boundary[side].edgegrid) ;

        if(!(Blocks[block].boundary[side].edgegrid 
        = (vertex *) malloc (Blocks[block].boundary[side].Imax*sizeof(vertex))))
    {
            printf("MEMORY ALLOCATION FAILURE 0 \n") ;
            return(1) ;
    }

    M = Blocks[block].boundary[side].Imax ;
    x1 = Blocks[block].corner[side].x ;
    x2 = Blocks[block].corner[NEXT(side)].x ;
    y1 = Blocks[block].corner[side].y ;
    y2 = Blocks[block].corner[NEXT(side)].y ;
    z1 = Blocks[block].corner[side].z ;
    z2 = Blocks[block].corner[NEXT(side)].z ;
    zi_stretching = Blocks[block].boundary[side].zi_stretching ;

        /* TFI interpolation for edge */
    for(i=0; i<M; i++)
    {
        zi = (i*1.0)/(M-1.0) ;

        if(side < 2)
        { 
            zi = pow(zi,zi_stretching) ;
        }
        else
        {
            zi = (1-pow((1-zi),zi_stretching)) ;
        }

        Blocks[block].boundary[side].edgegrid[i].x = (1-zi)*x1 + (zi)*x2 ;
        Blocks[block].boundary[side].edgegrid[i].y = (1-zi)*y1 + (zi)*y2 ;
        Blocks[block].boundary[side].edgegrid[i].z = (1-zi)*z1 + (zi)*z2 ;
    }
    Blocks[block].boundary[side].gridded = 1 ;

    return(0) ;
}

int generate_TFI_block(int block)
{
    long int    i, j,M,N ;
    double  zi, eta ;
    double  x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4 ;
    double  zi_stretching, eta_stretching ;

    if(Blocks[block].gridded == 1) return ;

//  free(Blocks[0].grid) ;

        if(!(Blocks[block].grid 
        = (vertex **) malloc (Blocks[block].Imax*sizeof(vertex *))))
    {
        printf("MEMORY ALLOCATION FAILURE 1 \n") ;
        return(1) ;
    }
 
        for(j=0; j<Blocks[block].Imax; j++)
    {
             if(!( Blocks[block].grid[j]
                          = (vertex *) malloc (Blocks[block].Jmax*sizeof(vertex))))
        {
            printf("j=%ld MEMORY ALLOCATION FAILURE 2 \n",j) ;
            return(2) ;
        }
    }

    M = Blocks[block].Imax ;
    N = Blocks[block].Jmax ;

    x1 = Blocks[block].boundary[0].edgegrid[0].x ;   
    x2 = Blocks[block].boundary[1].edgegrid[0].x ;   
    x3 = Blocks[block].boundary[2].edgegrid[0].x ;   
    x4 = Blocks[block].boundary[3].edgegrid[0].x ;   
    y1 = Blocks[block].boundary[0].edgegrid[0].y ;   
    y2 = Blocks[block].boundary[1].edgegrid[0].y ;   
    y3 = Blocks[block].boundary[2].edgegrid[0].y ;   
    y4 = Blocks[block].boundary[3].edgegrid[0].y ;   
    z1 = Blocks[block].boundary[0].edgegrid[0].z ;   
    z2 = Blocks[block].boundary[1].edgegrid[0].z ;   
    z3 = Blocks[block].boundary[2].edgegrid[0].z ;   
    z4 = Blocks[block].boundary[3].edgegrid[0].z ;   

    zi_stretching =  Blocks[block].zi_stretching ;
    eta_stretching =  Blocks[block].eta_stretching ;

        /*** TFI for 2D domain **********/
    for(i=0; i<M; i++)
    {
        zi = (i*1.0)/(M-1.0) ;

        zi = pow(zi,zi_stretching) ;

        for(j=0; j<N; j++)
        {
            eta = (j*1.0)/(N-1.0) ;

            eta = pow(eta,eta_stretching) ;
    
            Blocks[block].grid[i][j].x = 
                             (1-zi)*Blocks[block].boundary[3].edgegrid[N-1-j].x
                           +   (zi)*Blocks[block].boundary[1].edgegrid[j].x
                           + (1-eta)*Blocks[block].boundary[0].edgegrid[i].x
                           +   (eta)*Blocks[block].boundary[2].edgegrid[M-1-i].x
                           - (1-zi)*(1-eta)*x1 - (zi)*(1-eta)*x2
                           - (1-zi)*(eta)*x4   - (zi)*(eta)*x3 ;
            Blocks[block].grid[i][j].y = 
                             (1-zi)*Blocks[block].boundary[3].edgegrid[N-1-j].y
                           +   (zi)*Blocks[block].boundary[1].edgegrid[j].y
                           + (1-eta)*Blocks[block].boundary[0].edgegrid[i].y
                           +   (eta)*Blocks[block].boundary[2].edgegrid[M-1-i].y
                           - (1-zi)*(1-eta)*y1 - (zi)*(1-eta)*y2
                           - (1-zi)*(eta)*y4   - (zi)*(eta)*y3 ;
            Blocks[block].grid[i][j].z = 
                                                     (1-zi)*Blocks[block].boundary[3].edgegrid[N-1-j].z
                           +   (zi)*Blocks[block].boundary[1].edgegrid[j].z
                           + (1-eta)*Blocks[block].boundary[0].edgegrid[i].z
                           +   (eta)*Blocks[block].boundary[2].edgegrid[M-1-i].z
                           - (1-zi)*(1-eta)*z1 - (zi)*(1-eta)*z2
                           - (1-zi)*(eta)*z4   - (zi)*(eta)*z3 ;
        }
    }
    Blocks[block].gridded = 1 ;
    return(0) ;
}

void output_format()
{
    long int    i, j, ii, jj, M, N, count ;
    FILE        *fp ;
        long int        num, nCells, nVerts, nBounds ;
    long int    v1, v2, v3, v4 ;
    int         btype1, btype2, btype3 ;
    long int    nEdges1, nEdges2, nEdges3 ;

    fp = fopen("out.grd","w") ;

        nCells = 0 ;
        nVerts = 0 ;
    nBounds = 2 ;

        nCells += (Blocks[0].Imax-1) * (Blocks[0].Jmax-1) ; // block 1
        nVerts += (Blocks[0].Imax-1) * (Blocks[0].Jmax) ; // block 1
    
    fprintf(fp,"%d %d %d\n",nCells,nVerts,nBounds) ;
    printf("\n nCells = %d    nVerts = %d    nBounds = %d\n",nCells,nVerts,nBounds) ;

        /** writting out vertices **/

    /* block 1 */
    i = 0 ;
    for(jj=0; jj<Blocks[i].Jmax; jj++)
    for(ii=0; ii<Blocks[i].Imax-1; ii++)
      fprintf(fp,"%1.16e %1.16e\n",Blocks[i].grid[ii][jj].x,Blocks[i].grid[ii][jj].y) ;

        /* numbering vertices */
    count = 1 ;

    /* block 1 */
    i = 0 ;
    for(jj=0; jj<Blocks[i].Jmax; jj++)
    {
       for(ii=0; ii<Blocks[i].Imax-1; ii++)
       {
                Blocks[i].grid[ii][jj].vertexno = count++ ;
       }
           Blocks[i].grid[Blocks[i].Imax-1][jj].vertexno = Blocks[i].grid[0][jj].vertexno ;
    }

    /** writing cells **/

    for(i=0; i<Block_num; i++)
    {
    for(jj=0; jj<Blocks[i].Jmax-1; jj++)
    for(ii=0; ii<Blocks[i].Imax-1; ii++)
    {
        v1 = Blocks[i].grid[ii][jj].vertexno ;
        v2 = Blocks[i].grid[ii+1][jj].vertexno ;
        v3 = Blocks[i].grid[ii+1][jj+1].vertexno ;
        v4 = Blocks[i].grid[ii][jj+1].vertexno ;

        fprintf(fp,"%d %d %d %d\n",v1,v2,v3,v4) ;
    }
    }
    
    /** writing out boundaries **/      
    btype1 = 11 ;
    btype2 = 20 ;

    nEdges1 = nEdges2 = 0 ;

    nEdges1 += Blocks[0].Imax - 1 ;
    nEdges2 += Blocks[0].Imax - 1 ;

    /** Boundary 1 ** inner circle **/
    fprintf(fp,"%d %d\n",btype1,nEdges1) ;
    i = 0 ;jj = 0 ;
    for(ii=0; ii<Blocks[i].Imax-1; ii++)
    {
        v1 = Blocks[i].grid[ii][jj].vertexno ;
        v2 = Blocks[i].grid[ii+1][jj].vertexno ;
        fprintf(fp,"%d %d\n",v1,v2) ;
    }

    /** Boundary 2 **/
    fprintf(fp,"%d %d\n",btype2,nEdges2) ;
    i = 0 ; jj = Blocks[i].Jmax - 1 ;
    for(ii=Blocks[i].Imax-1; ii>0; ii--)
    {
        v1 = Blocks[i].grid[ii][jj].vertexno ;
        v2 = Blocks[i].grid[ii-1][jj].vertexno ;
        fprintf(fp,"%d %d\n",v1,v2) ;
    }
    
    fclose(fp) ;
    return ;
}

/*********************** output_plot3D *****************************************************/
void output_plot3D()
{       
        long int        i, j, ii, jj, M, N, count ;
        FILE            *fp ;
        long int        num, nCells, nVerts ;
        long int        v1, v2, v3, v4 ;
        int             izero ;
        double          dzero, dOne ;

        izero = 0 ;
        dzero = 0.00000000000000 ;
        dOne  = 1.00000000000000 ;


        fp = fopen("grid.x","w") ;
        
        nCells = 0 ;
        nVerts = 0 ;
        
        /* write number of blocks */
        fprintf(fp,"%d\n",1) ;
        printf("%d block(s)\n",1) ;
        
        nCells += (Blocks[0].Imax-1) * (Blocks[0].Jmax-1) ; // block 1
        nVerts += (Blocks[0].Imax) * (Blocks[0].Jmax) ; // block 1
                
        fprintf(fp,"%d %d %d\n",Blocks[0].Imax,Blocks[0].Jmax,2) ;
        printf("%d %d %d\n",Blocks[0].Imax,Blocks[0].Jmax,2) ;
        
        /** writting out vertices **/
        /* block 1 , x-coordinates */
        i = 0 ;
        for(jj=0; jj<Blocks[i].Jmax; jj++)
        for(ii=0; ii<Blocks[i].Imax; ii++)
        {
                fprintf(fp,"%1.16e\n",Blocks[i].grid[ii][jj].x) ;
        }
        for(jj=0; jj<Blocks[i].Jmax; jj++)
        for(ii=0; ii<Blocks[i].Imax; ii++)
        {
                fprintf(fp,"%1.16e\n",Blocks[i].grid[ii][jj].x) ;
        }
        /* block 1 , y-coordinates */
        i = 0 ;
        for(jj=0; jj<Blocks[i].Jmax; jj++)
        for(ii=0; ii<Blocks[i].Imax; ii++)
        {
                fprintf(fp,"%1.16e\n",Blocks[i].grid[ii][jj].y) ;
        }
        for(jj=0; jj<Blocks[i].Jmax; jj++)
        for(ii=0; ii<Blocks[i].Imax; ii++)
        {
                fprintf(fp,"%1.16e\n",Blocks[i].grid[ii][jj].y) ;
        }
        /* block 1 , z-coordinates */
        i = 0 ;
        for(jj=0; jj<Blocks[i].Jmax; jj++)
        for(ii=0; ii<Blocks[i].Imax; ii++)
        {
                fprintf(fp,"%1.16e\n",dzero) ;
        }
        for(jj=0; jj<Blocks[i].Jmax; jj++)
        for(ii=0; ii<Blocks[i].Imax; ii++)
        {
                fprintf(fp,"%1.16e\n",dOne) ;
        }

        fclose(fp) ;
        return ;
}

/*********************** output_tecplot ****************************************************/
void output_tecplot()
{
        long int        i, j, ii, jj, M, N, count ;
        FILE    *fp ;
        long int     num, nCells, nVerts, nBounds ;
        long int        v1, v2, v3, v4 ;
        int     btype1, btype2, btype3, btype4, btype5 ;
        long int        nEdges1, nEdges2, nEdges3, nEdges4, nEdges5 ;
        char    filename[15] ;

    nCells = 0 ;
        nVerts = 0 ;
        nBounds = 5 ;

        nCells += (Blocks[0].Imax-1) * (Blocks[0].Jmax-1) ; // block 1
        nVerts += (Blocks[0].Imax-1) * (Blocks[0].Jmax) ; // block 1

        /*** writing tecplot file ***/
        sprintf(filename,"out.plt") ;
        fp = fopen(filename,"w") ;

        /* writing header */
        fprintf(fp,"TITLE = \"GRID\"\n");
        fprintf(fp,"VARIABLES = \"X\",\"Y\"\n");
        fprintf(fp,"ZONE F = FEPOINT, N=%d, E=%d, ET=QUADRILATERAL\n",nVerts,nCells) ;

        /* writing node data */
    
    /* block 1 */
    i = 0 ;
    for(jj=0; jj<Blocks[i].Jmax; jj++)
    for(ii=0; ii<Blocks[i].Imax-1; ii++)
      fprintf(fp,"%1.16e %1.16e\n",Blocks[i].grid[ii][jj].x,Blocks[i].grid[ii][jj].y) ;

       /* numbering vertices */
    count = 1 ;

    /* block 1 */
    i = 0 ;
    for(jj=0; jj<Blocks[i].Jmax; jj++)
    {
       for(ii=0; ii<Blocks[i].Imax-1; ii++)
       {
                Blocks[i].grid[ii][jj].vertexno = count++ ;
       }
           Blocks[i].grid[Blocks[i].Imax-1][jj].vertexno = Blocks[i].grid[0][jj].vertexno ;
    }

        /* writing element data */

        for(i=0; i<Block_num; i++)
        {
        for(jj=0; jj<Blocks[i].Jmax-1; jj++)
        for(ii=0; ii<Blocks[i].Imax-1; ii++)
        {
                v1 = Blocks[i].grid[ii][jj].vertexno ;
                v2 = Blocks[i].grid[ii+1][jj].vertexno ;
                v3 = Blocks[i].grid[ii+1][jj+1].vertexno ;
                v4 = Blocks[i].grid[ii][jj+1].vertexno ;

                fprintf(fp,"%d %d %d %d\n",v1,v2,v3,v4) ;
        }
        }

        fclose(fp) ;
        /* wrote tecplot file */

}
/******************************************************************************************/

void testcase()
{
    vertex  A, B, C ;
    vertex  centre ;

    long int i, id, side, M ;
    int blk1, blk2, corner1, corner2 ;

    return ;
}
/******************************************************************************************/


/***************************************************************************************************/

int laplace_updatexyz(int iBlock, int i, int j)
{

    Blocks[1].grid[i][j].x = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].x 
                               + Blocks[iBlock].grid[i-1][j+1].x 
                               + Blocks[iBlock].grid[i+1][j-1].x 
                               + Blocks[iBlock].grid[i-1][j-1].x ) ; 
    
    Blocks[1].grid[i][j].y = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].y 
                               + Blocks[iBlock].grid[i-1][j+1].y 
                               + Blocks[iBlock].grid[i+1][j-1].y 
                               + Blocks[iBlock].grid[i-1][j-1].y ) ; 
    
    return 0 ;
}

int runlaplacian(int iBlock)
{
    int i, j, iter, iterate ;
    int error=0 ;
    
        iterate = 100 ;

        for(i=0;i<Imax;i++)
        {   
      for(j=0;j<Jmax;j++)
      {
        Blocks[1].grid[i][j].x = Blocks[iBlock].grid[i][j].x ; 
        Blocks[1].grid[i][j].y = Blocks[iBlock].grid[i][j].y ; 
      }
    }

    for(iter=0; iter<iterate; iter++)
    {
          for(i=1;i<Imax-1;i++)
          { 
        for(j=1;j<Jmax-1;j++)
        {
            Blocks[1].grid[i][j].x = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].x 
                               + Blocks[iBlock].grid[i-1][j+1].x 
                               + Blocks[iBlock].grid[i+1][j-1].x 
                               + Blocks[iBlock].grid[i-1][j-1].x ) ; 
    
            Blocks[1].grid[i][j].y = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].y 
                               + Blocks[iBlock].grid[i-1][j+1].y 
                               + Blocks[iBlock].grid[i+1][j-1].y 
                               + Blocks[iBlock].grid[i-1][j-1].y ) ; 
        }
          }

          i = 0 ;
        for(j=1;j<Jmax-1;j++)
        {
            Blocks[1].grid[i][j].x = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].x 
                               + Blocks[iBlock].grid[i+1][j+1].x 
                               + Blocks[iBlock].grid[i+1][j-1].x 
                               + Blocks[iBlock].grid[i+1][j-1].x ) ; 
    
            Blocks[1].grid[i][j].y = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].y 
                               + Blocks[iBlock].grid[i+1][j+1].y 
                               + Blocks[iBlock].grid[i+1][j-1].y 
                               + Blocks[iBlock].grid[i+1][j-1].y ) ; 
        }

          i = Imax-1 ;
        for(j=1;j<Jmax-1;j++)
        {
            Blocks[1].grid[i][j].x = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i-1][j+1].x 
                               + Blocks[iBlock].grid[i-1][j+1].x 
                               + Blocks[iBlock].grid[i-1][j-1].x 
                               + Blocks[iBlock].grid[i-1][j-1].x ) ; 
    
            Blocks[1].grid[i][j].y = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i-1][j+1].y 
                               + Blocks[iBlock].grid[i-1][j+1].y 
                               + Blocks[iBlock].grid[i-1][j-1].y 
                               + Blocks[iBlock].grid[i-1][j-1].y ) ; 
        }
          
          for(i=0;i<Imax;i++)
          { 
        for(j=0;j<Jmax;j++)
        {
          Blocks[iBlock].grid[i][j].x = Blocks[1].grid[i][j].x ; 
          Blocks[iBlock].grid[i][j].y = Blocks[1].grid[i][j].y ; 
        }
      }

        }
    
        return 0 ;
} 

int runlaplacian1(int iBlock)
{
    int i, j, iter, iterate ;
    int error=0 ;
    
        iterate = 100 ;

    for(iter=0; iter<iterate; iter++)
    {
          i = 0 ;
        for(j=1;j<Jmax-1;j++)
        {
            Blocks[0].grid[i][j].x = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].x 
                               + Blocks[iBlock].grid[i+1][j+1].x 
                               + Blocks[iBlock].grid[i+1][j-1].x 
                               + Blocks[iBlock].grid[i+1][j-1].x ) ; 
    
            Blocks[0].grid[i][j].y = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].y 
                               + Blocks[iBlock].grid[i+1][j+1].y 
                               + Blocks[iBlock].grid[i+1][j-1].y 
                               + Blocks[iBlock].grid[i+1][j-1].y ) ; 
        }

          for(i=1;i<Imax-1;i++)
          { 
        for(j=1;j<Jmax-1;j++)
        {
            Blocks[0].grid[i][j].x = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].x 
                               + Blocks[iBlock].grid[i-1][j+1].x 
                               + Blocks[iBlock].grid[i+1][j-1].x 
                               + Blocks[iBlock].grid[i-1][j-1].x ) ; 
    
            Blocks[0].grid[i][j].y = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i+1][j+1].y 
                               + Blocks[iBlock].grid[i-1][j+1].y 
                               + Blocks[iBlock].grid[i+1][j-1].y 
                               + Blocks[iBlock].grid[i-1][j-1].y ) ; 
        }
          }

          i = Imax-1 ;
        for(j=1;j<Jmax-1;j++)
        {
            Blocks[0].grid[i][j].x = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i-1][j+1].x 
                               + Blocks[iBlock].grid[i-1][j+1].x 
                               + Blocks[iBlock].grid[i-1][j-1].x 
                               + Blocks[iBlock].grid[i-1][j-1].x ) ; 
    
            Blocks[0].grid[i][j].y = (1.0/(4.0))*( 
                             Blocks[iBlock].grid[i-1][j+1].y 
                               + Blocks[iBlock].grid[i-1][j+1].y 
                               + Blocks[iBlock].grid[i-1][j-1].y 
                               + Blocks[iBlock].grid[i-1][j-1].y ) ; 
        }
        }
    
        return 0 ;
} 


