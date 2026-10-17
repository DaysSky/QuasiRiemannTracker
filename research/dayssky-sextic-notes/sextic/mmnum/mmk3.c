// Kernel: for each row u = x + y w and each spec s, compute sum_{i in [lo_s,hi_s)} w_s[i] chi_{n_i}(u).
// Sextic symbols via per-prime lookup tables of the exponent e (chi = zeta6^e), 6 = zero.
// Build: gcc -O3 -march=native -fopenmp mmk.c -o mmk -lm
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

static int64_t powmod(int64_t a, int64_t e, int64_t m){ __int128 r=1, b=a%m; if(b<0)b+=m; while(e){ if(e&1) r=(r*b)%m; b=(b*b)%m; e>>=1;} return (int64_t)r; }
typedef struct { int64_t a,b; } fq2;
static fq2 fmul(fq2 x, fq2 y, int64_t q){ fq2 r; r.a=((x.a*y.a - x.b*y.b)%q+q)%q; r.b=((x.a*y.b + x.b*y.a - x.b*y.b)%q+q)%q; return r; }
static fq2 fpow(fq2 x, int64_t e, int64_t q){ fq2 r={1,0}; while(e){ if(e&1) r=fmul(r,x,q); x=fmul(x,x,q); e>>=1;} return r; }

static void rd(void *p, size_t sz, size_t n, FILE *f){ if(fread(p,sz,n,f)!=n){ fprintf(stderr,"read error\n"); exit(1);} }

int main(int argc, char **argv){
  if(argc<3){ fprintf(stderr,"usage: mmk in.bin out.bin\n"); return 1; }
  FILE *f=fopen(argv[1],"rb"); int64_t hdr[4]; rd(hdr,8,4,f);
  int64_t np=hdr[0], ni=hdr[1], nr=hdr[2], ns=hdr[3];
  int64_t *kind=malloc(8*np), *P=malloc(8*np), *par=malloc(8*np);
  rd(kind,8,np,f); rd(P,8,np,f); rd(par,8,np,f);
  int64_t *spf=malloc(8*ni), *cof=malloc(8*ni); rd(spf,8,ni,f); rd(cof,8,ni,f);
  int64_t *rx=malloc(8*nr), *ry=malloc(8*nr); rd(rx,8,nr,f); rd(ry,8,nr,f);
  int64_t *lo=malloc(8*ns), *hi=malloc(8*ns); rd(lo,8,ns,f); rd(hi,8,ns,f);
  int64_t tw=0, *woff=malloc(8*ns); for(int64_t s=0;s<ns;s++){ woff[s]=tw; tw+=hi[s]-lo[s]; }
  double *w=malloc(8*(tw>0?tw:1)); rd(w,8,tw,f); fclose(f);
  int64_t maxhi=0; for(int64_t s=0;s<ns;s++) if(hi[s]>maxhi) maxhi=hi[s];
  // no tables: precompute for split primes the residue r and zeta6 image; inert: q
  int64_t *z6s=malloc(8*np*6);
  for(int64_t j=0;j<np;j++){
    if(kind[j]==0){ int64_t p=P[j], r=par[j]; int64_t z=1; for(int k=0;k<6;k++){ z6s[6*j+k]=z; z=(int64_t)(((__int128)z*((1+r)%p))%p); } }
  }
  uint8_t comb[49]; for(int a=0;a<7;a++) for(int b=0;b<7;b++) comb[a*7+b]=(a==6||b==6)?6:(a+b)%6;
  double cs[7], sn[7]; for(int k=0;k<6;k++){ cs[k]=cos(M_PI*k/3); sn[k]=sin(M_PI*k/3);} cs[6]=sn[6]=0;
  double *out=malloc(sizeof(double)*3*nr*ns);
  // primes needed: those with index < some bound; we just do all
  #pragma omp parallel
  {
    uint8_t *ep=malloc(np), *e=malloc(maxhi>0?maxhi:1);
    #pragma omp for schedule(dynamic,64)
    for(int64_t ir=0; ir<nr; ir++){
      int64_t x=rx[ir], y=ry[ir];
      for(int64_t j=0;j<np;j++){
        int64_t p=P[j];
        if(kind[j]==0){ int64_t t=(x+y*par[j])%p; if(t<0)t+=p; if(t==0){ep[j]=6;} else { int64_t h=powmod(t,(p-1)/6,p); int e=6; for(int k=0;k<6;k++) if(z6s[6*j+k]==h){e=k;break;} ep[j]=(uint8_t)e; } }
        else { int64_t q=par[j]; int64_t tx=x%q; if(tx<0)tx+=q; int64_t ty=y%q; if(ty<0)ty+=q;
               if(tx==0&&ty==0){ep[j]=6;} else { fq2 u={tx,ty}; fq2 h=fpow(u,(q*q-1)/6,q); fq2 zz={1,0}; fq2 z6={1,1}; int e=6; for(int k=0;k<6;k++){ if(zz.a==h.a&&zz.b==h.b){e=k;break;} zz=fmul(zz,z6,q);} ep[j]=(uint8_t)e; } }
      }
      e[0]=0; for(int64_t i=1;i<maxhi;i++) e[i]=comb[ep[spf[i]]*7+e[cof[i]]];
      for(int64_t s=0;s<ns;s++){
        double acc[7]={0,0,0,0,0,0,0}; const double *ws=w+woff[s]-lo[s];
        for(int64_t i=lo[s];i<hi[s];i++) acc[e[i]]+=ws[i];
        double re=0,im=0; for(int k=0;k<6;k++){ re+=acc[k]*cs[k]; im+=acc[k]*sn[k]; }
        double ms=0; for(int k=0;k<6;k++) ms+=acc[k]; out[3*(ir*ns+s)]=re; out[3*(ir*ns+s)+1]=im; out[3*(ir*ns+s)+2]=ms;
      }
    }
    free(ep); free(e);
  }
  FILE *g=fopen(argv[2],"wb"); fwrite(out,sizeof(double),3*nr*ns,g); fclose(g);
  return 0;
}
