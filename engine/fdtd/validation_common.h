#ifndef RFSYN_VALIDATION_COMMON_H
#define RFSYN_VALIDATION_COMMON_H
#ifndef RFSYN_VALIDATION_REFINEMENT
#define RFSYN_VALIDATION_REFINEMENT 1
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"
#include "materials.h"
static double pulse(double t, double t0, double tau) { double u=(t-t0)/tau; return -2*u*exp(-u*u); }
static double db(double x) { return x>0 ? 20*log10(x) : -INFINITY; }
static double sample(const fdtd_grid_t *g,size_t i) { return g->ez[fdtd_index(g,i,0,0)]; }
static double hsample(const fdtd_grid_t *g,size_t i) { return 0.5*(g->hy[fdtd_index(g,i,0,0)]+g->hy[fdtd_index(g,i-1,0,0)]); }
static void sheet(fdtd_grid_t *g,size_t i,double v) {
 for(size_t j=0;j<g->ny;j++) for(size_t k=0;k<g->nz;k++) g->ez[fdtd_index(g,i,j,k)]+=v;
}
/* modes: 0 CPML, 1 Fresnel, 2 energy. Both runs use production kernels. */
static int plane_validation_eps(int mode, double eps_r) {
 const size_t scale=RFSYN_VALIDATION_REFINEMENT;
 const size_t steps=1600*scale, src=100*scale, probe=180*scale, interface=240*scale, trans=260*scale;
 const size_t nx=(mode?800:320)*scale, npml=24*scale;
 int status=1; FILE *out=NULL; double *a=NULL,*b=NULL;
 fdtd_grid_t *r=NULL,*g=NULL; cpml_t *pr=NULL,*pg=NULL; materials_t *m=NULL;
 r=fdtd_grid_create(1000*scale,4,4,1e-3/scale,0.99); g=fdtd_grid_create(nx,4,4,1e-3/scale,0.99);
 if(!r||!g) goto cleanup;
 if(fdtd_grid_set_periodic(r,6)||fdtd_grid_set_periodic(g,6)) goto cleanup;
 pr=cpml_create(r,cpml_default_params(npml)); pg=cpml_create(g,cpml_default_params(npml));
 m=materials_create(g); a=calloc(steps,sizeof(*a)); b=calloc(steps,sizeof(*b));
 if(!pr||!pg||!m||!a||!b) goto cleanup;
 if(mode && materials_set_box(m,g,interface,nx+1,0,5,0,5,eps_r,0)) goto cleanup;
 /* Interface lies on the integer-x E plane: tangential dual cell is half each medium. */
 if(mode) for(size_t j=0;j<4;j++) for(size_t k=0;k<4;k++) {
  if(materials_set_component(m,g,FIELD_EZ,interface,j,k,0.5*(1+eps_r),0) || materials_set_component(m,g,FIELD_EY,interface,j,k,0.5*(1+eps_r),0)) goto cleanup;
 }
 const char *name=mode==0?"cpml_planewave_validation.csv":mode==1?"materials_fresnel_validation.csv":(eps_r==1?"materials_energy_vacuum_control.csv":"materials_energy_validation.csv");
 out=fopen(name,"w"); if(!out) goto cleanup;
 if(mode==2) fprintf(out,"step,t,energy_incident_cum,energy_reflected_cum,energy_transmitted_cum\n");
 double inc=0,refl=0,tr=0,peak=0,err=0,ss=0,rr=0;
 for(size_t n=0;n<steps;n++) {
  double hr=hsample(r,probe), hg=hsample(g,probe), ht=hsample(g,trans);
  fdtd_update_h_cpml(r,pr); fdtd_update_h_cpml(g,pg);
  double er=sample(r,probe), eg=sample(g,probe), et=sample(g,trans);
  inc+=-er*0.5*(hr+hsample(r,probe))*g->dt;
  refl+=-(eg-er)*0.5*(hg+hsample(g,probe)-hr-hsample(r,probe))*g->dt;
  tr+=-et*0.5*(ht+hsample(g,trans))*g->dt;
  if(mode==2 && fprintf(out,"%zu,%.17g,%.17g,%.17g,%.17g\n",n,n*g->dt,inc,refl,tr)<0) goto cleanup;
  fdtd_update_e_cpml(r,pr); fdtd_update_e_cpml_materials(g,pg,m);
  /* Injection after E update: same vacuum location and timing in both runs. */
  double v=pulse((n+1)*g->dt,60*scale*g->dt,16*scale*g->dt); sheet(r,src,v); sheet(g,src,v);
  a[n]=sample(r,probe); b[n]=sample(g,probe);
  if(!isfinite(a[n])||!isfinite(b[n])||!isfinite(inc)||!isfinite(refl)||!isfinite(tr)) goto cleanup;
  if(fabs(a[n])>peak) peak=fabs(a[n]);
 }
 if(!(peak>1e-12)||!(inc>0)) goto cleanup;
 if(mode!=2) fprintf(out,"step,t,probe_ref,probe_test,abs_err,db_err\n");
 double reflected_peak=0;
 for(size_t n=0;n<steps;n++) {
  double d=fabs(b[n]-a[n]); if(d>err) err=d; ss+=d*d; rr+=a[n]*a[n];
  if(mode && n>=350*scale && n<=600*scale && d>reflected_peak) reflected_peak=d;
  if(mode!=2 && fprintf(out,"%zu,%.17g,%.17g,%.17g,%.17g,%.17g\n",n,(n+1)*g->dt,a[n],b[n],d,db(d/peak))<0) goto cleanup;
 }
 if(mode==0) {
  printf("CPML: max normalized error %.9g (%.3f dB), relative L2 %.9g\n",err/peak,db(err/peak),sqrt(ss/rr));
  status=err/peak<=1e-3?0:1;
 } else if(mode==1) {
  double theory=fabs((1-sqrt(eps_r))/(1+sqrt(eps_r)));
  double measured=reflected_peak/peak, rel=fabs(measured-theory)/theory;
  printf("Fresnel pulse magnitude: measured %.9g, theory %.9g, relative error %.6g%%\n",measured,theory,100*rel);
  status=rel<=0.02?0:1;
 } else {
  double balance=(-refl+tr)/inc, mismatch=fabs(balance-1);
  printf("Energy per unit area: R=%.9g T=%.9g balance=%.9g mismatch=%.6g%%\n",-refl/inc,tr/inc,balance,100*mismatch);
  status=(refl<=0 && tr>=0 && mismatch<=0.02 && (eps_r!=1 || fabs(refl/inc)<=1e-12))?0:1;
 }
cleanup:
 if(out && fclose(out)) status=1;
 free(a);free(b);materials_destroy(m);cpml_destroy(pr);cpml_destroy(pg);fdtd_grid_destroy(r);fdtd_grid_destroy(g);
 printf("%s\n",status?"FAIL":"PASS"); return status;
}
static int plane_validation(int mode) {
 if(mode==2 && plane_validation_eps(2,1)) return 1;
 return plane_validation_eps(mode,4);
}
#endif
