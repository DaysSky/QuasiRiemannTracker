// gen.c  -- enumerate squarefree primary c in Z[i], N(c) <= Z, and compute
//   G(c)  = g(conj chi_c) = sum_{x mod c} conj(chi_c(x)) psi(x/c),  psi(z)=exp(2 pi i Im z),
//   e_k(c) = exponent of chi_c(m_k) (chi_c = (./c)_4) for the test primes m_k (255 = zero),
// via twisted multiplicativity G(c'pi) = conj chi_{c'}(pi) conj chi_pi(c') G(c') G(pi)
//   = s * (c'/pi)_2 * G(c') G(pi),  s = (-1)^{((Nc'-1)/4)((Npi-1)/4)}  (quartic reciprocity).
// usage: gen Z NTEST out.bin      (test primes = all primary primes of norm <= NTEST)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>
typedef long long ll;
typedef unsigned char u8;
#define MAXK 64
typedef struct { ll a,b,N,p,r; int inert; double complex G; u8 e[MAXK]; } prime_t;

static ll mulmod(ll a, ll b, ll m){ return (ll)((__int128)a*b % m); }
static ll powmod(ll b, ll e, ll m){ ll r=1; b%=m; if(b<0)b+=m; while(e){ if(e&1) r=mulmod(r,b,m); b=mulmod(b,b,m); e>>=1;} return r; }
static ll inv(ll a, ll m){ return powmod(a,m-2,m); }
static int is_primary(ll a, ll b){ ll A=((a%4)+4)%4, B=((b%4)+4)%4; return (A==1&&B==0)||(A==3&&B==2); }
static void prim_assoc(ll a, ll b, ll *pa, ll *pb){
  for(int u=0;u<4;u++){ if(is_primary(a,b)){*pa=a;*pb=b;return;} ll t=a; a=-b; b=t; }
  fprintf(stderr,"no primary assoc\n"); exit(1);
}
static ll primroot(ll p){
  ll f[64]; int nf=0; ll m=p-1;
  for(ll q=2;q*q<=m;q++) if(m%q==0){ f[nf++]=q; while(m%q==0) m/=q; }
  if(m>1) f[nf++]=m;
  for(ll g=2;;g++){ int ok=1; for(int i=0;i<nf;i++) if(powmod(g,(p-1)/f[i],p)==1){ok=0;break;} if(ok) return g; }
}
// Gaussian integer mod q (inert) power
static void gpowq(ll a, ll b, ll e, ll q, ll *ra, ll *rb){
  ll x=1,y=0; a%=q; if(a<0)a+=q; b%=q; if(b<0)b+=q;
  while(e){ if(e&1){ ll nx=(mulmod(x,a,q)-mulmod(y,b,q)+q)%q, ny=(mulmod(x,b,q)+mulmod(y,a,q))%q; x=nx;y=ny; }
    ll na=(mulmod(a,a,q)-mulmod(b,b,q)+q)%q, nb=mulmod(2*a%q,b,q); a=na;b=nb; e>>=1; }
  *ra=x;*rb=y;
}
static int unit_exp_q(ll x, ll y, ll q){ // which unit i^k is (x,y) mod q
  if(x==1&&y==0) return 0; if(x==0&&y==1) return 1; if(x==q-1&&y==0) return 2; if(x==0&&y==q-1) return 3; return -1;
}

static prime_t *P; static int np=0;
static ll Zmax; static int K;
// output
static ll *oN,*oa,*ob; static double *oGr,*oGi; static u8 *oe; static ll no=0, cap=0;
static void emit(ll a, ll b, ll N, double complex G, const u8 *e){
  if(no>=cap){ cap=cap?2*cap:1<<20; oN=realloc(oN,cap*8); oa=realloc(oa,cap*8); ob=realloc(ob,cap*8);
    oGr=realloc(oGr,cap*8); oGi=realloc(oGi,cap*8); oe=realloc(oe,cap*K); }
  oN[no]=N; oa[no]=a; ob[no]=b; oGr[no]=creal(G); oGi[no]=cimag(G); memcpy(oe+no*K,e,K); no++;
}
static void dfs(int start, ll a, ll b, ll N, double complex G, const u8 *e){
  emit(a,b,N,G,e);
  u8 e2[MAXK];
  for(int t=start;t<np;t++){
    prime_t *q=&P[t];
    if(N*q->N>Zmax) break;
    int leg;
    if(!q->inert){ ll y=((a%q->p)+q->p + mulmod(((b%q->p)+q->p)%q->p,q->r,q->p))%q->p; ll z=powmod(y,(q->p-1)/2,q->p); leg=(z==1)?1:-1; if(y==0){fprintf(stderr,"y=0\n");exit(1);} }
    else { ll y=N% q->p; ll z=powmod(y,(q->p-1)/2,q->p); leg=(z==1)?1:-1; if(y==0){fprintf(stderr,"y=0 inert\n");exit(1);} }
    int s=((((N-1)/4)&1)&&(((q->N-1)/4)&1))?-1:1;
    double complex G2=(double)(s*leg)*G*q->G;
    for(int k=0;k<K;k++) e2[k]=(e[k]==255||q->e[k]==255)?255:((e[k]+q->e[k])&3);
    ll na=a*q->a-b*q->b, nb=a*q->b+b*q->a;
    dfs(t+1,na,nb,N*q->N,G2,e2);
  }
}
static int cmpp(const void *x,const void *y){ const prime_t *A=x,*B=y; if(A->N!=B->N) return A->N<B->N?-1:1; if(A->a!=B->a) return A->a<B->a?-1:1; return A->b<B->b?-1:(A->b>B->b); }

int main(int argc,char**argv){
  Zmax=atoll(argv[1]); ll NT=atoll(argv[2]); const char*out=argv[3];
  // sieve
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
  K=0; while(K<np && P[K].N<=NT) K++;
  if(K>MAXK){fprintf(stderr,"K too big\n");return 1;}
  fprintf(stderr,"np=%d K=%d\n",np,K);
  // per-prime data
  #pragma omp parallel
  {
    u8 *tab=NULL; ll tabcap=0;
    #pragma omp for schedule(dynamic,16)
    for(int t=0;t<np;t++){
      prime_t *q=&P[t];
      if(!q->inert){
        ll p=q->p, A=((q->a%p)+p)%p, B=((q->b%p)+p)%p;
        ll r=mulmod(p-A,inv(B,p),p); q->r=r;
        ll g=primroot(p); ll z=powmod(g,(p-1)/4,p); int ee = (z==r)?1:((z==p-r)?3:-1);
        if(ee<0){fprintf(stderr,"bad e p=%lld\n",p); exit(1);}
        if(p>tabcap){ tab=realloc(tab,p); tabcap=p; }
        ll y=1; for(ll s=0;s<p-1;s++){ tab[y]=(u8)((ee*s)&3); y=mulmod(y,g,p); }
        // G = sum_x i^{-k(x)} exp(-2 pi i B x/p)
        double Sr[4]={0,0,0,0}, Si[4]={0,0,0,0};
        double ang=-2*M_PI*(double)B/(double)p; double complex w=cexp(I*ang), zz=1.0;
        for(ll x=1;x<p;x++){
          if((x&1023)==0){ ll ph=mulmod(B,x,p); zz=cexp(-2*M_PI*I*(double)ph/(double)p); } else zz*=w;
          int k=tab[x]; Sr[k]+=creal(zz); Si[k]+=cimag(zz);
        }
        // i^{-k}: k=0:1, k=1:-i, k=2:-1, k=3:i
        double complex S0=Sr[0]+I*Si[0],S1=Sr[1]+I*Si[1],S2=Sr[2]+I*Si[2],S3=Sr[3]+I*Si[3];
        q->G = S0 - I*S1 - S2 + I*S3;
        for(int k=0;k<K;k++){
          ll ma=P[k].a, mb=P[k].b;
          ll yy=(((ma%p)+p)%p + mulmod(((mb%p)+p)%p,r,p))%p;
          if(yy==0){ q->e[k]=255; continue; }
          ll zz2=powmod(yy,(p-1)/4,p);
          q->e[k]= (zz2==1)?0:(zz2==r)?1:(zz2==p-1)?2:(zz2==p-r)?3:254;
          if(q->e[k]==254){fprintf(stderr,"bad chi\n");exit(1);}
        }
      } else {
        ll qq=q->p; double complex G=0; ll E=(qq*qq-1)/4;
        for(ll u=0;u<qq;u++) for(ll v=0;v<qq;v++){ if(!u&&!v) continue; ll x,y; gpowq(u,v,E,qq,&x,&y); int k=unit_exp_q(x,y,qq);
          // conj chi = i^{-k};  psi(x/(-q)) = exp(-2 pi i v/q)   [pi=-q]
          double complex c = (k==0)?1:(k==1)?-I:(k==2)?-1:I;
          G += c*cexp(-2*M_PI*I*(double)v/(double)qq); }
        // note: primary associate of q is -q (q = 3 mod 4); formula above assumes pi=-q
        if(!(q->a==-qq && q->b==0)){fprintf(stderr,"inert assoc unexpected\n");exit(1);}
        q->G=G;
        for(int k=0;k<K;k++){ ll ma=P[k].a, mb=P[k].b; ll x0=((ma%qq)+qq)%qq, y0=((mb%qq)+qq)%qq;
          if(x0==0&&y0==0){ q->e[k]=255; continue; }
          ll x,y; gpowq(x0,y0,E,qq,&x,&y); int kk=unit_exp_q(x,y,qq); if(kk<0){fprintf(stderr,"bad unit\n");exit(1);} q->e[k]=kk; }
      }
    }
    free(tab);
  }
  fprintf(stderr,"gauss sums done\n");
  u8 e0[MAXK]; memset(e0,0,MAXK);
  dfs(0,1,0,1,1.0,e0);
  fprintf(stderr,"nodes=%lld\n",no);
  FILE*f=fopen(out,"wb"); int32_t hdr[2]={(int32_t)K,0}; fwrite(hdr,4,2,f); fwrite(&no,8,1,f);
  fwrite(oN,8,no,f); fwrite(oa,8,no,f); fwrite(ob,8,no,f); fwrite(oGr,8,no,f); fwrite(oGi,8,no,f); fwrite(oe,1,no*K,f);
  // test primes info
  for(int k=0;k<K;k++){ ll v[3]={P[k].a,P[k].b,P[k].N}; fwrite(v,8,3,f); double g[2]={creal(P[k].G),cimag(P[k].G)}; fwrite(g,8,2,f); }
  fclose(f);
  return 0;
}
