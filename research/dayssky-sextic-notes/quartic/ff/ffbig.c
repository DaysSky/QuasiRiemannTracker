// Function-field (L3) in the large-degree regime: q prime = 1 mod 4, rows m monic deg d, columns c monic squarefree deg n.
// T(c) = chi_c(-1) * (c'/c)_4 (up to a unit depending only on n), R(m) = sum_c T(c) conj chi_c(m).
// Symbols via factorisation of c into irreducibles: (a/pi)_4 = a^{(q^k-1)/4} mod pi in F_q, mapped to i^t.
// usage: ffbig q n d
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#define MAXD 16
typedef struct { int deg; int c[MAXD+1]; } poly; // monic unless noted, coeffs low->high
static int q, g0, dlogt[1024], qm14;
static int md(long x){ x%=q; if(x<0) x+=q; return (int)x; }
static void pmodp(poly *a, const poly *b){ // a := a mod b (b monic)
  while(a->deg>=b->deg && a->deg>=0){ int f=a->c[a->deg]; if(f){ int s=a->deg-b->deg; for(int i=0;i<=b->deg;i++) a->c[s+i]=md(a->c[s+i]-(long)f*b->c[i]); }
    while(a->deg>=0 && a->c[a->deg]==0) a->deg--; }
}
static void pmulmod(poly *r, const poly *a, const poly *b, const poly *m){
  poly t; t.deg=a->deg+b->deg; if(a->deg<0||b->deg<0){ r->deg=-1; return; }
  memset(t.c,0,sizeof(t.c)); long tmp[2*MAXD+2]; memset(tmp,0,sizeof(tmp));
  for(int i=0;i<=a->deg;i++) for(int j=0;j<=b->deg;j++) tmp[i+j]+= (long)a->c[i]*b->c[j];
  // reduce tmp mod m
  int deg=t.deg; for(int i=0;i<=deg;i++) tmp[i]%=q;
  while(deg>=m->deg){ long f=tmp[deg]%q; if(f){ int s=deg-m->deg; for(int i=0;i<=m->deg;i++) tmp[s+i]=(tmp[s+i]-f*m->c[i])%q; } deg--; }
  r->deg=-1; for(int i=0;i<=deg && i<m->deg;i++){ r->c[i]=md(tmp[i]); if(r->c[i]) r->deg=i; }
  for(int i=r->deg+1;i<=MAXD;i++) r->c[i]=0;
}
static int quart(const poly *a0, const poly *pi){ // (a/pi)_4 as exponent t in Z/4 (i^t), or -1 if zero
  poly a=*a0; pmodp(&a,pi); if(a.deg<0) return -1;
  long e=1; for(int i=0;i<pi->deg;i++) e*=q; e=(e-1)/4;
  poly r; r.deg=0; memset(r.c,0,sizeof(r.c)); r.c[0]=1; poly b=a;
  while(e){ if(e&1) pmulmod(&r,&r,&b,pi); pmulmod(&b,&b,&b,pi); e>>=1; }
  if(r.deg!=0){ fprintf(stderr,"bad quartic\n"); exit(1); }
  return dlogt[r.c[0]]/qm14;
}
// enumerate monic irreducibles up to degree n via sieve on all monic polys (encode as base-q integers)
static long enc(const poly *a){ long v=0; for(int i=a->deg-1;i>=0;i--) v=v*q+a->c[i]; return v; }
static void dec(long v,int deg,poly *a){ a->deg=deg; for(int i=0;i<deg;i++){ a->c[i]=v%q; v/=q; } a->c[deg]=1; for(int i=deg+1;i<=MAXD;i++) a->c[i]=0; }
static void pmul(poly *r,const poly*a,const poly*b){ poly t; memset(&t,0,sizeof t); t.deg=a->deg+b->deg; for(int i=0;i<=a->deg;i++) for(int j=0;j<=b->deg;j++) t.c[i+j]=md(t.c[i+j]+(long)a->c[i]*b->c[j]); *r=t; }
int main(int argc,char**argv){
  q=atoi(argv[1]); int n=atoi(argv[2]), d=atoi(argv[3]);
  for(g0=2;g0<q;g0++){ int ok=1; long x=1; for(int k=1;k<q-1;k++){ x=x*g0%q; if(x==1){ok=0;break;} } if(ok) break; }
  { long x=1; for(int k=0;k<q-1;k++){ dlogt[x]=k; x=x*g0%q; } } qm14=(q-1)/4;
  int N=n>d?n:d;
  // factor table: for each monic poly of degree <= N store smallest irreducible factor index via sieve
  long *cnt=malloc((N+1)*sizeof(long)); long tot=0; for(int k=0;k<=N;k++){ long v=1; for(int i=0;i<k;i++) v*=q; cnt[k]=v; }
  // irreducibles of degree k: polys with no factor of degree <= k/2; brute-force sieve: mark products
  char **isred=malloc((N+1)*sizeof(char*)); for(int k=1;k<=N;k++){ isred[k]=calloc(cnt[k],1); }
  poly *irr=NULL; int nirr=0, cap=0;
  for(int k=1;k<=N;k++){
    for(long v=0;v<cnt[k];v++) if(!isred[k][v]){ if(nirr>=cap){cap=cap?2*cap:1024; irr=realloc(irr,cap*sizeof(poly));} dec(v,k,&irr[nirr]); nirr++; }
    // mark products of irreducible of degree k with all monic of degree j (k+j<=N)
    for(int idx=0;idx<nirr;idx++) if(irr[idx].deg==k){
      for(int j=1;j+k<=N;j++) for(long w=0;w<cnt[j];w++){ poly b,pr; dec(w,j,&b); pmul(&pr,&irr[idx],&b); isred[k+j][enc(&pr)/1]=1; }
    }
  }
  fprintf(stderr,"irreducibles up to degree %d: %d\n",N,nirr);
  // factorisation by trial division using irreducibles (deg <= n)
  // columns: monic squarefree c of degree n with factor list
  typedef struct { int nf; int f[MAXD]; double complex T; } col;
  col *C=malloc(cnt[n]*sizeof(col)); long nc=0;
  poly minus1; minus1.deg=0; memset(minus1.c,0,sizeof minus1.c); minus1.c[0]=q-1;
  #pragma omp parallel for schedule(dynamic,256)
  for(long v=0;v<cnt[n];v++){
    poly c; dec(v,n,&c); poly rem=c; int nf=0, f[MAXD]; int sq=1;
    for(int idx=0;idx<nirr && rem.deg>0;idx++){ if(irr[idx].deg>rem.deg) continue;
      poly t=rem; pmodp(&t,&irr[idx]); if(t.deg>=0) continue;
      // divide
      int times=0; while(1){ poly t2=rem; pmodp(&t2,&irr[idx]); if(t2.deg>=0) break;
        // exact division rem/irr
        poly qq; memset(&qq,0,sizeof qq); qq.deg=rem.deg-irr[idx].deg; poly a=rem;
        for(int s=qq.deg;s>=0;s--){ int co=a.c[s+irr[idx].deg]; qq.c[s]=co; for(int i=0;i<=irr[idx].deg;i++) a.c[s+i]=md(a.c[s+i]-(long)co*irr[idx].c[i]); }
        rem=qq; times++; }
      if(times>1) sq=0; f[nf++]=idx; }
    if(!sq) { C[v].nf=-1; continue; }
    // T(c) = chi_c(-1) (c'/c)_4
    poly der; memset(&der,0,sizeof der); der.deg=n-1; for(int i=1;i<=n;i++) der.c[i-1]=md((long)i*c.c[i]); while(der.deg>=0&&der.c[der.deg]==0) der.deg--;
    int ex=0, zero=0; for(int k=0;k<nf;k++){ int a=quart(&minus1,&irr[f[k]]); int b=quart(&der,&irr[f[k]]); if(a<0||b<0){zero=1;break;} ex+=a+b; }
    C[v].nf=nf; for(int k=0;k<nf;k++) C[v].f[k]=f[k];
    ex&=3; C[v].T = zero?0:((ex==0)?1:(ex==1)?I:(ex==2)?-1:-I);
  }
  // rows m: all monic of degree d; chi_c(m) = prod over c's factors pi of (m/pi)_4 -> need table m x irreducibles(deg<=n)
  long NM=cnt[d];
  // table[m][idx] exponent (char, -1 if zero)
  signed char *tab=malloc((size_t)NM*nirr);
  #pragma omp parallel for schedule(dynamic,16)
  for(long w=0;w<NM;w++){ poly m; dec(w,d,&m); for(int idx=0;idx<nirr;idx++){ if(irr[idx].deg>n){ tab[(size_t)w*nirr+idx]=0; continue;} int e=quart(&m,&irr[idx]); tab[(size_t)w*nirr+idx]=(signed char)e; } }
  double S=0, D=0, S0=0; long rows=0; double maxr=0; double S_sqf=0,D_sqf=0;
  #pragma omp parallel for reduction(+:S,D,S0,S_sqf,D_sqf) reduction(max:maxr) schedule(dynamic,16)
  for(long w=0;w<NM;w++){
    double complex R=0, R0=0; double dg=0;
    for(long v=0;v<cnt[n];v++){ if(C[v].nf<0) continue; int ex=0, zero=0;
      for(int k=0;k<C[v].nf;k++){ int e=tab[(size_t)w*nirr+C[v].f[k]]; if(e<0){zero=1;break;} ex+=e; }
      if(zero) continue; ex&=3; double complex z=(ex==0)?1:(ex==1)?-I:(ex==2)?-1:I; // conj chi
      R+=C[v].T*z; R0+=z; dg+=1; }
    S+=creal(R)*creal(R)+cimag(R)*cimag(R); D+=dg; double rr=(creal(R)*creal(R)+cimag(R)*cimag(R))/(dg>0?dg:1); if(rr>maxr) maxr=rr;
    // squarefree rows only
    poly m; dec(w,d,&m); poly der; memset(&der,0,sizeof der); der.deg=d-1; for(int i=1;i<=d;i++) der.c[i-1]=md((long)i*m.c[i]); while(der.deg>=0&&der.c[der.deg]==0) der.deg--;
    // gcd(m,m') == 1 ?
    poly a=m,b=der; while(b.deg>=0){ // make b monic copy for pmodp
      int inv=1; for(int t=1;t<q;t++) if((long)t*b.c[b.deg]%q==1){inv=t;break;} poly bm=b; for(int i=0;i<=bm.deg;i++) bm.c[i]=md((long)bm.c[i]*inv); poly t2=a; pmodp(&t2,&bm); a=bm; b=t2; }
    if(a.deg==0){ S_sqf+=creal(R)*creal(R)+cimag(R)*cimag(R); D_sqf+=dg; }
  }
  printf("q=%d n=%d d=%d rows=%ld  S/diag=%.4f  squarefree-rows S/diag=%.4f  max row ratio=%.2f\n",q,n,d,NM,S/D,S_sqf/D_sqf,maxr);
  return 0;
}
