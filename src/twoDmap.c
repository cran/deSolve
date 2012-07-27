/* --------------------------------------------------------------------*
 SPARSITY of 2-D and 3-D reaction-transport problems with mapping
 the states that are present have a value > 0 in vector 'ipres' 
 ipres contains the actual number of state variable, 
 after applying the mask , e.g. ipres(20) = 10 means that the element
 20 in the original 2D matrix is the 10th element, after applying the mask
  -------------------------------------------------------------------- */
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include <R_ext/Rdynload.h>
#include "deSolve.h"


void sparsity2Dmap (SEXP Type, int* iwork, int neq, int liw) {

    int nspec, nx, ny, bndx, bndy, Nt, ij, isp, i, j, k, l, m;
    int totN, *ipres, Mnew;
    nspec = INTEGER(Type)[1]; /* number components*/
    nx    = INTEGER(Type)[2]; /* dimension x*/
    ny    = INTEGER(Type)[3]; /* dimension y*/
    bndx  = INTEGER(Type)[4]; /* cyclic boundary x*/
    bndy  = INTEGER(Type)[5]; /* cyclic boundary y*/
    totN  = INTEGER(Type)[7]; /* Total state variables in original 2D matrix*/
    
	  ipres = (int *) R_alloc(totN, sizeof(int));
    for (j=0; j < totN; j++) ipres[j] = INTEGER(Type)[j+8];
	
    Nt    = nx*ny;
    ij    = 31 + neq;
    iwork[30] = 1;
    m = 1;
    for( i = 0; i < nspec; i++) {
      isp = i*Nt;
      for( j = 0; j < nx; j++) {
        for( k = 0; k < ny; k++) {
          if (ij > liw-8-nspec)  
            error("not enough memory allocated in iwork - increase liw %i ",liw);
          Mnew = ipres[m-1];
          if (Mnew > 0) {
          
             interactmap (&ij, liw, iwork, ipres, m);
             if (k < ny-1)     interactmap (&ij, liw, iwork, ipres, m+1); 
             if (j < nx-1)     interactmap (&ij, liw, iwork, ipres, m+ny);
             if (j > 0)        interactmap (&ij, liw, iwork, ipres, m-ny);
             if (k > 0)        interactmap (&ij, liw, iwork, ipres, m-1);
             if (bndx == 1) {
               if (j == 0)     interactmap (&ij, liw, iwork, ipres, isp+(nx-1)*ny+k+1);
               if (j == nx-1)  interactmap (&ij, liw, iwork, ipres, isp+k+1);
             }
             if (bndy == 1) {
               if (k == 0)    interactmap (&ij, liw, iwork, ipres, isp+(j+1)*ny);
               if (k == ny-1) interactmap (&ij, liw, iwork, ipres, isp + j*ny +1);
             }
             for(l = 0; l < nspec; l++)
               if (l != i)    interactmap (&ij, liw, iwork, ipres, l*Nt+j*ny+k+1);

             iwork[30+Mnew] = ij-30-neq;
          }
          m = m+1;
        }
      }
    }
}

void interactmap (int *ij, int nnz, int *iwork, int *ipres, int ival) {
  
/* check if not yet present for current state */
     if (ipres[ival-1] > 0) {
       if (*ij > nnz) 
         error ("not enough memory allocated in iwork - increase liw %i ", nnz);
     iwork[(*ij)++] = ipres[ival-1];

  }
}


/*==================================================*/
/* an element in C-array A(I,J,K), i=0,dim(1)-1 etc... is positioned at 
   j*dim(2)*dim(3) + k*dim(3) + l + 1 in FORTRAN VECTOR! 
   includes check on validity

   dimens and boundary are reversed ... 
*/

void sparsity3Dmap (SEXP Type, int* iwork, int neq, int liw) {
    int nspec, nx, ny, nz, bndx, bndy, bndz, Nt, ij, isp, i, j, k, l, m, ll;
    int totN, *ipres, Mnew;

    nspec = INTEGER(Type)[1]; 
    nx    = INTEGER(Type)[2]; 
    ny    = INTEGER(Type)[3]; 
    nz    = INTEGER(Type)[4]; 
    bndx  = INTEGER(Type)[5];
    bndy  = INTEGER(Type)[6]; 
    bndz  = INTEGER(Type)[7]; 
    totN  = INTEGER(Type)[9]; /* Total state variables in original 3D matrix*/

	ipres = (int *) R_alloc(totN, sizeof(int));
     for (j=0; j < totN; j++) {ipres[j] = INTEGER(Type)[j+10];
     }

    Nt    = nx*ny*nz;
    ij    = 31+neq;
    iwork[30] = 1;
    m = 1;
    for( i = 0; i < nspec; i++) {
      isp = i*Nt;
      for( j = 0; j < nx; j++) {
        for( k = 0; k < ny; k++) {
          for( ll = 0; ll < nz; ll++) {
            
            if (ij > liw-6-nspec)  
              error ("not enough memory allocated in iwork - increase liw %i ", liw);

            Mnew = ipres[m-1];
            if (Mnew > 0) {
              interactmap (&ij, liw, iwork, ipres, m);
              if (ll < nz-1)  
                interactmap (&ij, liw, iwork, ipres, m+1);
              else if (bndz == 1)
                interactmap (&ij, liw, iwork, ipres, isp + j*ny*nz + k*nz + 1);              
            
              if (k  < ny-1) 
                interactmap (&ij, liw, iwork, ipres, m+nz);
              else  if (bndy == 1)
                interactmap (&ij, liw, iwork, ipres, isp + j*ny*nz + ll + 1);
               
              if (j  < nx-1)  
                interactmap (&ij, liw, iwork, ipres, m+ny*nz);
              else if (bndx == 1)  
                interactmap (&ij, liw, iwork, ipres, isp + k*nz + ll + 1);

              if (j > 0)      
                interactmap (&ij, liw, iwork, ipres, m-ny*nz);
              else if (bndx == 1)
                interactmap (&ij, liw, iwork, ipres, isp+(nx-1)*ny*nz+k*nz+ll+1);
                                         
              if (k > 0)  
                interactmap (&ij, liw, iwork, ipres, m-nz);
              else  if (bndy == 1)
                interactmap (&ij, liw, iwork, ipres, isp + j*ny*nz+(ny-1)*nz+ll+1);      
            
              if (ll > 0) 
                interactmap (&ij, liw, iwork, ipres, m-1);
              else if (bndz == 1)  
                interactmap (&ij, liw, iwork, ipres, isp + j*ny*nz+k*nz+nz);

              for(l = 0; l < nspec; l++)
                if (l != i) 
                  interactmap (&ij, liw, iwork, ipres, l*Nt+j*ny*nz+k*nz+ll+1);
              iwork[30+Mnew] = ij-30-neq;
            }
            m = m+1;
          }
        }
      }
    }
}

