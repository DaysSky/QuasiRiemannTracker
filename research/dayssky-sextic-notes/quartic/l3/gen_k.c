// gen_rows.c -- rows of the (G) mean square for the quartic family over Z[i].
//   R(m) = sum_{c sqfree primary, Z/2 < N(c) <= Z} T(c) conj(chi_c(m)) W(N(c)/Z),
//   T(c) = alpha(c) * gtilde(conj chi_c)   (twisted row; mode 1)
//   T0(c) = gtilde(conj chi_c)             (untwisted control; mode 0)
// for all primary m with N(m) <= Q.  Gauss sums and quartic symbols as in gen_big.c
// (twisted multiplicativity G(c' pi) = s (c'/pi)_2 G(c') G(pi), s = (-1)^{((Nc'-1)/4)((Npi-1)/4)}).
// Output: summary line + per-row |R|^2 (both modes) as binary.
// usage: gen_rows Z Q out.bin
#include "helpers.inc"
#include <omp.h>

typedef struct { ll a,b,N,p,r; int inert; double complex G; } prime_t;
static prime_t *P; static int np=0; static u8 *PE; // PE[t*K+k] = exponent of chi_{P_t}(m_k)
static ll Zmax, Qmax; static int K;
// rows
static int NR; static int *rfac_n; static int (*rfac)[24]; static ll *rN;
// batch of c
#define BATCH 65536
static double complex *bw1,*bw0; static u8 *be; static int nb=0;
static double *S1,*S0;            // per-row accumulators reduced at end
static double complex *R1,*R0;    // global rows (complex)
static double diagw=0, *rowdiag;  // diagonal sum_c |W|^2 over c coprime to m (per row)
static double wbump(double x){ if(x<=0.5||x>=1.0) return 0; return exp(-1.0/((x-0.5)*(1.0-x))+16.0); }

static void flush(void){
  if(!nb) return;
  #pragma omp parallel
  {
    double complex *l1=calloc(NR,sizeof(double complex)), *l0=calloc(NR,sizeof(double complex));
    double *ld=calloc(NR,sizeof(double));
    #pragma omp for schedule(static)
    for(int j=0;j<nb;j++){
      const u8 *e=be+(ll)j*K; double complex w1=bw1[j], w0=bw0[j]; double ww=creal(w0)*creal(w0)+cimag(w0)*cimag(w0);
      for(int m=0;m<NR;m++){
        int ex=0, zero=0;
        for(int f=0;f<rfac_n[m];f++){ int k=rfac[m][f]; u8 v=e[k]; if(v==255){zero=1;break;} ex+=v; }
        if(zero) continue;
        ex&=3; // conj chi = i^{-ex}
        double complex z = (ex==0)?1.0:(ex==1)?I:(ex==2)?-1.0:-I;
        l1[m]+=w1*z; l0[m]+=w0*z; ld[m]+=ww;
      }
    }
    #pragma omp critical
    { for(int m=0;m<NR;m++){ R1[m]+=l1[m]; R0[m]+=l0[m]; rowdiag[m]+=ld[m]; } }
    free(l1); free(l0); free(ld);
  }
  nb=0;
}
static void emitc(ll a, ll b, ll N, double complex G, const u8 *e, int mu){
  double x=(double)N/(double)Zmax; double W=wbump(x); if(W==0) return;
  double complex gt=G/sqrt((double)N);
  double complex al=((double)a+I*(double)b)/sqrt((double)N);
  bw1[nb]=mu*W*gt; bw0[nb]=W*gt; (void)al; memcpy(be+(ll)nb*K,e,K); nb++;
  if(nb==BATCH) flush();
}
static void dfs(int start, ll a, ll b, ll N, double complex G, const u8 *e, int mu){
  if(2*N>Zmax) emitc(a,b,N,G,e,mu);
  u8 *e2=malloc(K);
  for(int t=start;t<np;t++){
    prime_t *q=&P[t];
    if(N*q->N>Zmax) break;
    int leg;
    if(!q->inert){ ll y=((a%q->p)+q->p + mulmod(((b%q->p)+q->p)%q->p,q->r,q->p))%q->p; ll z=powmod(y,(q->p-1)/2,q->p); leg=(z==1)?1:-1; }
    else { ll y=N% q->p; ll z=powmod(y,(q->p-1)/2,q->p); leg=(z==1)?1:-1; }
    int s=((((N-1)/4)&1)&&(((q->N-1)/4)&1))?-1:1;
    double complex G2=(double)(s*leg)*G*q->G;
    const u8 *qe=PE+(ll)t*K;
    for(int k=0;k<K;k++) e2[k]=(e[k]==255||qe[k]==255)?255:((e[k]+qe[k])&3);
    ll na=a*q->a-b*q->b, nb2=a*q->b+b*q->a;
    dfs(t+1,na,nb2,N*q->N,G2,e2,-mu);
  }
  free(e2);
}
static int cmpp(const void *x,const void *y){ const prime_t *A=x,*B=y; if(A->N!=B->N) return A->N<B->N?-1:1; if(A->a!=B->a) return A->a<B->a?-1:1; return A->b<B->b?-1:(A->b>B->b); }
// rows: all primary m (products of test primes with multiplicity), N(m) <= Q; m=1 included
static int rcap=0;
static void rowdfs(int start, ll N, int *fl, int nf){
  if(NR>=rcap){ rcap=rcap?2*rcap:4096; rfac_n=realloc(rfac_n,rcap*sizeof(int)); rfac=realloc(rfac,rcap*sizeof(*rfac)); rN=realloc(rN,rcap*8); }
  rfac_n[NR]=nf; for(int i=0;i<nf;i++) rfac[NR][i]=fl[i]; rN[NR]=N; NR++;
  for(int k=start;k<K;k++){ if(N*P[k].N>Qmax) break; if(nf>=24) break; fl[nf]=k; rowdfs(k,N*P[k].N,fl,nf+1); }
}
int main(int argc,char**argv){
  Zmax=atoll(argv[1]); Qmax=atoll(argv[2]); const char*out=argv[3];
  char *isc=calloc(Zmax+1,1); for(ll i=2;i*i<=Zmax;i++) if(!isc[i]) for(ll j=i*i;j<=Zmax;j+=i) isc[j]=1;
  ll cnt=0; for(ll p=3;p<=Zmax;p++) if(!isc[p]){ if(p%4==1) cnt+=2; else if(p*p<=Zmax) cnt++; }
  P=calloc(cnt,sizeof(prime_t));
  for(ll p=3;p<=Zmax;p++) if(!isc[p]){
    if(p%4==1){ ll a=0,b=0; for(a=1;a*a<p;a++){ ll b2=p-a*a; ll s=(ll)llround(sqrt((double)b2)); if(s*s==b2){b=s;break;} }
      ll A,B; prim_assoc(a,b,&A,&B); P[np].a=A;P[np].b=B;P[np].N=p;P[np].p=p;P[np].inert=0;np++;
      P[np].a=A;P[np].b=-B;P[np].N=p;P[np].p=p;P[np].inert=0;np++; }
    else if(p*p<=Zmax){ ll A,B; prim_assoc(p,0,&A,&B); P[np].a=A;P[np].b=B;P[np].N=p*p;P[np].p=p;P[np].inert=1;np++; }
  }
  qsort(P,np,sizeof(prime_t),cmpp);
  K=0; while(K<np && P[K].N<=Qmax) K++;
  PE=malloc((ll)np*K);
  fprintf(stderr,"np=%d K=%d\n",np,K);
  #pragma omp parallel
  {
    u8 *tab=NULL; ll tabcap=0;
    #pragma omp for schedule(dynamic,16)
    for(int t=0;t<np;t++){
      prime_t *q=&P[t]; u8 *qe=PE+(ll)t*K;
      if(!q->inert){
        ll p=q->p, A=((q->a%p)+p)%p, B=((q->b%p)+p)%p;
        ll r=mulmod(p-A,inv(B,p),p); q->r=r;
        ll g=primroot(p); ll z=powmod(g,(p-1)/4,p); int ee = (z==r)?1:((z==p-r)?3:-1);
        if(ee<0){fprintf(stderr,"bad e p=%lld\n",p); exit(1);}
        if(p>tabcap){ tab=realloc(tab,p); tabcap=p; }
        ll y=1; for(ll s=0;s<p-1;s++){ tab[y]=(u8)((ee*s)&3); y=mulmod(y,g,p); }
        double Sr[4]={0,0,0,0}, Si[4]={0,0,0,0};
        double ang=-2*M_PI*(double)B/(double)p; double complex w=cexp(I*ang), zz=1.0;
        for(ll x=1;x<p;x++){
          if((x&1023)==0){ ll ph=mulmod(B,x,p); zz=cexp(-2*M_PI*I*(double)ph/(double)p); } else zz*=w;
          int k=tab[x]; Sr[k]+=creal(zz); Si[k]+=cimag(zz);
        }
        double complex S0=Sr[0]+I*Si[0],S1v=Sr[1]+I*Si[1],S2=Sr[2]+I*Si[2],S3=Sr[3]+I*Si[3];
        q->G = S0 - I*S1v - S2 + I*S3;
        for(int k=0;k<K;k++){
          ll ma=P[k].a, mb=P[k].b;
          ll yy=(((ma%p)+p)%p + mulmod(((mb%p)+p)%p,r,p))%p;
          if(yy==0){ qe[k]=255; continue; }
          ll zz2=powmod(yy,(p-1)/4,p);
          qe[k]= (zz2==1)?0:(zz2==r)?1:(zz2==p-1)?2:(zz2==p-r)?3:254;
          if(qe[k]==254){fprintf(stderr,"bad chi\n");exit(1);}
        }
      } else {
        ll qq=q->p; double complex G=0; ll E=(qq*qq-1)/4;
        for(ll u=0;u<qq;u++) for(ll v=0;v<qq;v++){ if(!u&&!v) continue; ll x,y; gpowq(u,v,E,qq,&x,&y); int k=unit_exp_q(x,y,qq);
          double complex c = (k==0)?1:(k==1)?-I:(k==2)?-1:I;
          G += c*cexp(-2*M_PI*I*(double)v/(double)qq); }
        if(!(q->a==-qq && q->b==0)){fprintf(stderr,"inert assoc unexpected\n");exit(1);}
        q->G=G;
        for(int k=0;k<K;k++){ ll ma=P[k].a, mb=P[k].b; ll x0=((ma%qq)+qq)%qq, y0=((mb%qq)+qq)%qq;
          if(x0==0&&y0==0){ qe[k]=255; continue; }
          ll x,y; gpowq(x0,y0,E,qq,&x,&y); int kk=unit_exp_q(x,y,qq); if(kk<0){fprintf(stderr,"bad unit\n");exit(1);} qe[k]=kk; }
      }
    }
    free(tab);
  }
  int fl[24]; rowdfs(0,1,fl,0);
  fprintf(stderr,"gauss sums done; rows NR=%d\n",NR);
  R1=calloc(NR,sizeof(double complex)); R0=calloc(NR,sizeof(double complex)); rowdiag=calloc(NR,sizeof(double));
  bw1=malloc(BATCH*sizeof(double complex)); bw0=malloc(BATCH*sizeof(double complex)); be=malloc((ll)BATCH*K);
  u8 *e0=calloc(K,1);
  dfs(0,1,0,1,1.0,e0,1); flush();
  double S1t=0,S0t=0,Dt=0;
  FILE*f=fopen(out,"wb"); int32_t hdr[2]={NR,0}; fwrite(hdr,4,2,f);
  for(int m=0;m<NR;m++){ double a1=creal(R1[m])*creal(R1[m])+cimag(R1[m])*cimag(R1[m]); double a0=creal(R0[m])*creal(R0[m])+cimag(R0[m])*cimag(R0[m]);
    S1t+=a1; S0t+=a0; Dt+=rowdiag[m];
    double rec[4]={(double)rN[m],a1,a0,rowdiag[m]}; fwrite(rec,8,4,f); }
  fclose(f);
  printf("{\"Z\":%lld,\"Q\":%lld,\"rows\":%d,\"K\":%d,\"S_twisted\":%.6e,\"S_untwisted\":%.6e,\"diag\":%.6e,\"ratio_twisted\":%.5f,\"ratio_untwisted\":%.5f}\n",
         Zmax,Qmax,NR,K,S1t,S0t,Dt,S1t/Dt,S0t/Dt);
  return 0;
}
