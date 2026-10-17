// primes_g X out.bin : for split primary primes pi (N=p<=X, p=1 mod 4), G = g(conj chi_pi) (same convention as gen_rows.c)
#include "helpers.inc"
#include <omp.h>
int main(int argc,char**argv){
  ll X=atoll(argv[1]); char*isc=calloc(X+1,1); for(ll i=2;i*i<=X;i++) if(!isc[i]) for(ll j=i*i;j<=X;j+=i) isc[j]=1;
  int n=0; ll *A=malloc(8*X/2),*B=malloc(8*X/2),*Pp=malloc(8*X/2);
  for(ll p=5;p<=X;p++) if(!isc[p]&&p%4==1){ ll a,b=0; for(a=1;a*a<p;a++){ ll b2=p-a*a; ll s=llround(sqrt((double)b2)); if(s*s==b2){b=s;break;} }
    ll aa,bb; prim_assoc(a,b,&aa,&bb); A[n]=aa;B[n]=bb;Pp[n]=p;n++; A[n]=aa;B[n]=-bb;Pp[n]=p;n++; }
  double *out=malloc(sizeof(double)*5*n);
  #pragma omp parallel
  { u8*tab=NULL; ll cap=0;
  #pragma omp for schedule(dynamic,8)
  for(int t=0;t<n;t++){ ll p=Pp[t], a=((A[t]%p)+p)%p, b=((B[t]%p)+p)%p; ll r=mulmod(p-a,inv(b,p),p);
    ll g=primroot(p); ll z=powmod(g,(p-1)/4,p); int ee=(z==r)?1:((z==p-r)?3:-1); if(ee<0){fprintf(stderr,"bad\n");exit(1);}
    if(p>cap){tab=realloc(tab,p);cap=p;} ll y=1; for(ll s=0;s<p-1;s++){ tab[y]=(u8)((ee*s)&3); y=mulmod(y,g,p);} 
    double Sr[4]={0},Si[4]={0}; double ang=-2*M_PI*(double)b/(double)p; double complex w=cexp(I*ang),zz=1.0;
    for(ll x=1;x<p;x++){ if((x&1023)==0){ ll ph=mulmod(b,x,p); zz=cexp(-2*M_PI*I*(double)ph/(double)p);} else zz*=w; int k=tab[x]; Sr[k]+=creal(zz); Si[k]+=cimag(zz);} 
    double complex G=(Sr[0]+I*Si[0]) - I*(Sr[1]+I*Si[1]) - (Sr[2]+I*Si[2]) + I*(Sr[3]+I*Si[3]);
    out[5*t]=A[t]; out[5*t+1]=B[t]; out[5*t+2]=p; out[5*t+3]=creal(G); out[5*t+4]=cimag(G); }
  free(tab);}
  FILE*f=fopen(argv[2],"wb"); fwrite(out,8,5*(ll)n,f); fclose(f); fprintf(stderr,"n=%d\n",n); return 0; }
