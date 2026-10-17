// ms.c -- mean square  sum_{N(m) in bin} |B_m(Z)|^2,  B_m = sum_c a_c conj(chi_c(m)) ,  c sqfree primary, Z/2<N(c)<=Z
// coefficient sets: 0: a_c = alpha(c) G~(c) W   (the (G) rows)
//                   1: a_c = G~(c) W            (untwisted)
//                   2,3: a_c = random phase * W
// rows: all nonzero Gaussian integers m with N(m) <= Qmax, binned by N(m) in [Qmax/2^(b+1), Qmax/2^b)
// usage: ms Z Qmax nbins seed
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
typedef long long ll; typedef unsigned char u8;
typedef struct { ll a,b,N,p,r; int inert; double complex G; u8 *tab; } prime_t;
static ll mulmod(ll a, ll b, ll m){ return (ll)((__int128)a*b % m); }
static ll powmod(ll b, ll e, ll m){ ll r=1; b%=m; if(b<0)b+=m; while(e){ if(e&1) r=mulmod(r,b,m); b=mulmod(b,b,m); e>>=1;} return r; }
static int is_primary(ll a, ll b){ ll A=((a%4)+4)%4, B=((b%4)+4)%4; return (A==1&&B==0)||(A==3&&B==2); }
static void prim_assoc(ll a, ll b, ll *pa, ll *pb){ for(int u=0;u<4;u++){ if(is_primary(a,b)){*pa=a;*pb=b;return;} ll t=a; a=-b; b=t; } exit(1); }
static ll primroot(ll p){ ll f[64]; int nf=0; ll m=p-1; for(ll q=2;q*q<=m;q++) if(m%q==0){ f[nf++]=q; while(m%q==0) m/=q; } if(m>1) f[nf++]=m;
  for(ll g=2;;g++){ int ok=1; for(int i=0;i<nf;i++) if(powmod(g,(p-1)/f[i],p)==1){ok=0;break;} if(ok) return g; } }
static void gpowq(ll a, ll b, ll e, ll q, ll *ra, ll *rb){ ll x=1,y=0; a%=q; if(a<0)a+=q; b%=q; if(b<0)b+=q;
  while(e){ if(e&1){ ll nx=(x*a-y*b)%q, ny=(x*b+y*a)%q; x=(nx+q)%q;y=(ny+q)%q; } ll na=(a*a-b*b)%q, nb=(2*a*b)%q; a=(na+q)%q;b=(nb+q)%q; e>>=1; } *ra=x;*rb=y; }
static int unit_exp_q(ll x, ll y, ll q){ if(x==1&&y==0) return 0; if(x==0&&y==1) return 1; if(x==q-1&&y==0) return 2; if(x==0&&y==q-1) return 3; return -1; }
static prime_t *P; static int np=0; static ll Z;
#define MAXW 12
typedef struct { ll a,b,N; double complex G; int nf; int f[MAXW]; } col_t;
static col_t *C; static ll nc=0, ccap=0;
static void dfs(int start, ll a, ll b, ll N, double complex G, int nf, int *f){
  if(2*N>Z){ if(nc>=ccap){ccap=ccap?2*ccap:4096; C=realloc(C,ccap*sizeof(col_t));} col_t *c=&C[nc++]; c->a=a;c->b=b;c->N=N;c->G=G;c->nf=nf; memcpy(c->f,f,nf*sizeof(int)); }
  for(int t=start;t<np;t++){ prime_t *q=&P[t]; if(N*q->N>Z) break; int leg;
    if(!q->inert){ ll y=((a%q->p)+q->p + mulmod(((b%q->p)+q->p)%q->p,q->r,q->p))%q->p; leg=(powmod(y,(q->p-1)/2,q->p)==1)?1:-1; }
    else { leg=(powmod(N%q->p,(q->p-1)/2,q->p)==1)?1:-1; }
    int s=((((N-1)/4)&1)&&(((q->N-1)/4)&1))?-1:1;
    f[nf]=t; dfs(t+1,a*q->a-b*q->b,a*q->b+b*q->a,N*q->N,(double)(s*leg)*G*q->G,nf+1,f); }
}
static int cmpp(const void *x,const void *y){ const prime_t *A=x,*B=y; if(A->N!=B->N) return A->N<B->N?-1:1; if(A->a!=B->a) return A->a<B->a?-1:1; return A->b<B->b?-1:(A->b>B->b); }
static double Wf(double x){ if(x<=0.5||x>=1) return 0; double y=4*x-3; return exp(1-1/(1-y*y)); }
int main(int argc,char**argv){
  Z=atoll(argv[1]); ll Qmax=atoll(argv[2]); int nb=atoi(argv[3]); unsigned seed=atoi(argv[4]);
  char *isc=calloc(Z+1,1); for(ll i=2;i*i<=Z;i++) if(!isc[i]) for(ll j=i*i;j<=Z;j+=i) isc[j]=1;
  ll cnt=0; for(ll p=3;p<=Z;p++) if(!isc[p]){ if(p%4==1) cnt+=2; else if(p*p<=Z) cnt++; }
  P=calloc(cnt,sizeof(prime_t));
  for(ll p=3;p<=Z;p++) if(!isc[p]){
    if(p%4==1){ ll a,b=0; for(a=1;a*a<p;a++){ ll b2=p-a*a; ll s=(ll)llround(sqrt((double)b2)); if(s*s==b2){b=s;break;} }
      ll A,B; prim_assoc(a,b,&A,&B); P[np]=(prime_t){A,B,p,p,0,0,0,NULL}; np++; P[np]=(prime_t){A,-B,p,p,0,0,0,NULL}; np++; }
    else if(p*p<=Z){ ll A,B; prim_assoc(p,0,&A,&B); P[np]=(prime_t){A,B,p*p,p,0,1,0,NULL}; np++; } }
  qsort(P,np,sizeof(prime_t),cmpp);
  #pragma omp parallel for schedule(dynamic,8)
  for(int t=0;t<np;t++){ prime_t *q=&P[t];
    if(!q->inert){ ll p=q->p, A=((q->a%p)+p)%p, B=((q->b%p)+p)%p; ll r=mulmod(p-A,powmod(B,p-2,p),p); q->r=r;
      ll g=primroot(p); ll z=powmod(g,(p-1)/4,p); int ee=(z==r)?1:3; q->tab=malloc(p); q->tab[0]=255;
      ll y=1; for(ll s=0;s<p-1;s++){ q->tab[y]=(u8)((ee*s)&3); y=mulmod(y,g,p); }
      double complex S[4]={0,0,0,0};
      for(ll x=1;x<p;x++){ S[q->tab[x]]+=cexp(-2*M_PI*I*(double)mulmod(B,x,p)/(double)p); }
      q->G=S[0]-I*S[1]-S[2]+I*S[3];
    } else { ll qq=q->p; ll E=(qq*qq-1)/4; q->tab=malloc(qq*qq); double complex G=0;
      for(ll u=0;u<qq;u++) for(ll v=0;v<qq;v++){ if(!u&&!v){ q->tab[0]=255; continue;} ll x,y; gpowq(u,v,E,qq,&x,&y); int k=unit_exp_q(x,y,qq); q->tab[u*qq+v]=k;
        double complex c=(k==0)?1:(k==1)?-I:(k==2)?-1:I; G+=c*cexp(-2*M_PI*I*(double)v/(double)qq); }
      q->G=G; } }
  int f[MAXW]; dfs(0,1,0,1,1.0,0,f);
  fprintf(stderr,"np=%d nc=%lld\n",np,nc);
  int NS=4; double complex *A=malloc(sizeof(double complex)*NS*nc); srand(seed); double sumsq[4]={0,0,0,0};
  for(ll i=0;i<nc;i++){ double w=Wf((double)C[i].N/(double)Z); double sq=sqrt((double)C[i].N); double complex al=(C[i].a+I*C[i].b)/sq, Gt=C[i].G/sq;
    A[0*nc+i]=al*Gt*w; A[1*nc+i]=Gt*w; A[2*nc+i]=cexp(2*M_PI*I*(rand()/(double)RAND_MAX))*w; A[3*nc+i]=cexp(2*M_PI*I*(rand()/(double)RAND_MAX))*w;
    for(int s=0;s<NS;s++) sumsq[s]+=creal(A[s*nc+i]*conj(A[s*nc+i])); }
  // rows
  ll R=(ll)sqrt((double)Qmax)+1;
  double *acc=calloc(nb*NS,sizeof(double)); double *mx=calloc(nb*NS,sizeof(double)); ll *rc=calloc(nb,sizeof(ll));
  // also record per-row values for top-row inspection (alpha set)
  #pragma omp parallel
  { u8 *em=malloc(np); double *lacc=calloc(nb*NS,sizeof(double)); ll *lrc=calloc(nb,sizeof(ll)); double *lmx=calloc(nb*NS,sizeof(double));
    #pragma omp for schedule(dynamic,4)
    for(ll ma=-R;ma<=R;ma++) for(ll mb=-R;mb<=R;mb++){ ll Nm=ma*ma+mb*mb; if(Nm==0||Nm>Qmax) continue;
      int bin=0; double qq=(double)Qmax; while(bin<nb && !( (double)Nm>qq/2 )){ qq/=2; bin++; } if(bin>=nb) continue;
      for(int t=0;t<np;t++){ prime_t *q=&P[t];
        if(!q->inert){ ll p=q->p; ll y=(((ma%p)+p)%p + mulmod(((mb%p)+p)%p,q->r,p))%p; em[t]=q->tab[y]; }
        else { ll qv=q->p; em[t]=q->tab[(((ma%qv)+qv)%qv)*qv+(((mb%qv)+qv)%qv)]; } }
      double complex Bs[4]={0,0,0,0};
      static const double complex IP[4]={1,-I,-1,I}; // conj chi = i^{-e}
      for(ll i=0;i<nc;i++){ int e=0,z=0; for(int j=0;j<C[i].nf;j++){ u8 v=em[C[i].f[j]]; if(v==255){z=1;break;} e+=v; } if(z) continue;
        double complex ph=IP[e&3]; for(int s=0;s<NS;s++) Bs[s]+=A[s*nc+i]*ph; }
      for(int s=0;s<NS;s++){ double v=creal(Bs[s]*conj(Bs[s])); lacc[bin*NS+s]+=v; if(v>lmx[bin*NS+s]) lmx[bin*NS+s]=v; } lrc[bin]++; }
    #pragma omp critical
    { for(int k=0;k<nb*NS;k++){ acc[k]+=lacc[k]; if(lmx[k]>mx[k]) mx[k]=lmx[k]; } for(int k=0;k<nb;k++) rc[k]+=lrc[k]; }
  }
  printf("Z=%lld Qmax=%lld nc=%lld  sumsq: alpha %.1f untw %.1f rand %.1f %.1f\n",Z,Qmax,nc,sumsq[0],sumsq[1],sumsq[2],sumsq[3]);
  printf("  Q-bin (N(m) in (Q/2,Q])   logQ/logZ  #rows   alpha   untw   rand1  rand2   [mean |B|^2 / sum|a|^2]\n");
  double qq=(double)Qmax;
  for(int b=0;b<nb;b++){ printf("  %10.0f  %6.3f  %8lld  %7.3f %7.3f %7.3f %7.3f\n",qq,log(qq)/log((double)Z),rc[b],
      acc[b*NS+0]/(rc[b]*sumsq[0]),acc[b*NS+1]/(rc[b]*sumsq[1]),acc[b*NS+2]/(rc[b]*sumsq[2]),acc[b*NS+3]/(rc[b]*sumsq[3]));
    printf("        max|B|^2/sum|a|^2:  alpha %.2f untw %.2f rand %.2f %.2f\n",mx[b*NS]/sumsq[0],mx[b*NS+1]/sumsq[1],mx[b*NS+2]/sumsq[2],mx[b*NS+3]/sumsq[3]); qq/=2; }
  return 0;
}
