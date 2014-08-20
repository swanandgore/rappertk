/* Numerical Recipes Sort EigenValues */
#include <math.h>
#include <vector>
#include<cstring>
#include<cstdlib>
#include<stdio.h>

using namespace std;

void EigenSort(vector<float> & values, vector<vector<float> > & vectors) {
    int k,j,i;
    float p;
    for (i=0;i<3;i++) {
        p = values[k=i];
        for (j=i+1;j<=3;j++) { if (values[j] >= p) {p=values[k=j];}}
        if (k != i) {
            values[k]=values[i];
            values[i]=p;
            for (j=0;j<=3;j++) {
                p=vectors[j][i];
                vectors[j][i]=vectors[j][k];
                vectors[j][k]=p;
            }
        }
    }
}


bool Jacobi (vector<vector<float> > & q, vector<float> & values, vector<vector<float> > & vectors) {
#define JACOBI_NUM 4
#define JACOBI_SQRNUM 16
#define JACOBI_MAXITER 9000
#define JACOBI_CONVERGENCE 1.0E-15
    int ip,iq, i,j, nrot;
    float tresh,theta,tau,t,s,h,g,c, b[4],z[4];
    double sm;

    s=0;
    z[0]=0.0; z[1]=0.0; z[2]=0.0; z[3]=0.0;
    vectors[0][0]= vectors[1][1]= vectors[2][2]= vectors[3][3]=1.0;
    vectors[0][1]= vectors[0][2]= vectors[0][3]= vectors[1][0]= vectors[1][2]= vectors[1][3]= vectors[2][0]= vectors[2][1]= vectors[2][3]= vectors[3][0]= vectors[3][1]= vectors[3][2]= 0.0;
    b[0]=q[0][0]; b[1]=q[1][1]; b[2]=q[2][2]; b[3]=q[3][3];
    values[0]=q[0][0]; values[1]=q[1][1]; values[2]=q[2][2]; values[3]=q[3][3];

    //printf("iteration 0 in routine jacobi sm=undefined values[0]=%f values[1]=%f values[2]=%f values[3]=%f \n",values[0],values[1],values[2],values[3]);
    nrot=0;
    /*JACOBI_MAXITER is number of iterations to do before giving up */
    i=1;
    while (i < JACOBI_MAXITER ) {
        sm=0.0; ip=0;
        while (ip<JACOBI_NUM-1) {
            iq=ip+1;
            while (iq<JACOBI_NUM) { sm += fabs(q[ip][iq]);iq++; }
            ip++;
        }
        //if (options.verbose > 4) printf("iteration %d in routine jacobi sm=%.28f values[0]=%f values[1]=%f values[2]=%f values[3]=%f \n",i,sm,values[0],values[1],values[2],values[3]);
        if (isnan (sm) != 0)  { printf("Cannot superimpose\n");return true;} /* For VERY different ( and large) structures*/
        /*  Rather than if sm == 0.0  - allow almost zero due to compiler rounding errors */
        if ( sm < JACOBI_CONVERGENCE ) { EigenSort(values,vectors); return false ; }
        if (i < 4) { tresh=0.2*sm/(float)(JACOBI_SQRNUM); }
        else       { tresh=0.0;  }
        ip=0;
        while (ip<JACOBI_NUM-1) {
            iq=ip+1;
            while  (iq<JACOBI_NUM)  {
                g=fabs(q[ip][iq])*(int)100;
                if ( i > 4 &&  fabs(values[ip]) > g  &&  fabs(values[iq]) > g  )  q[ip][iq]=0.0;
                else if (fabs(q[ip][iq]) > tresh ) {
                    h=values[iq]-values[ip];
                    if ( fabs(h)>g   ) {t= (q[ip][iq])/h;}
                    else {
                        theta =0.5*h/(q[ip][iq]);
                        t=  1.0/(fabs(theta)+sqrt((float)( 1.0+theta*theta)))  ;
                        if (theta < 0.0) t=-t;
                    }
                    c=1.0/sqrt((float)(1+t*t));
                    s = t*c;
                    tau = s/(1.0+c);
                    h = t*q[ip][iq];
                    z[ip] -= h;
                    z[iq] += h;
                    q[ip][iq] = 0.0;
                    values[ip] -= h;
                    values[iq] += h;
                    /* Case of rotations 0<=j<p */
                    j=0;
                    while (j< ip) { 
                        g=q[j][ip];h=q[j][iq];q[j][ip]=g-s*(h+g*tau); q[j][iq]=h+s*(g-h*tau);
                        j++; 
                    }
                    /* Case of rotations p<j<q */
                    j=ip+1;
                    while (j< iq) { 
                        g=q[ip][j];h=q[j][iq];q[ip][j]=g-s*(h+g*tau); q[j][iq]=h+s*(g-h*tau);
                        j++; 
                    }
                    /* Case of rotations q<j<=n */
                    j=iq+1;
                    while (j<JACOBI_NUM) { 
                        g=q[ip][j];h=q[iq][j];q[ip][j]=g-s*(h+g*tau); q[iq][j]=h+s*(g-h*tau);
                        j++; 
                    }
                    j=0;
                    while (j<JACOBI_NUM) {
                        g=vectors[j][ip];h=vectors[j][iq];vectors[j][ip]=g-s*(h+g*tau); vectors[j][iq]=h+s*(g-h*tau);
                        j++; 
                    }
                    nrot++;
                } /* else if */
                iq++;
            } /* for iq */
            ip++;
        } /* for ip */
        ip=0;
        while (ip<JACOBI_NUM) {
            b[ip]+=z[ip];
            values[ip]=b[ip];
            z[ip]=0.0;
            ip++;
        }
        i++;
    }  /* end of i loop */
    fprintf(stderr,"Too many iterations in routine jacobi sm=%f\n",sm);
#undef JACOBI_ROTATE
#undef JACOBI_NUM
#undef JACOBI_SQRNUM
#undef JACOBI_MAXITER
    return true;
}
