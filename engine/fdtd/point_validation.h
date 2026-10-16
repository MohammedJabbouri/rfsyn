#ifndef RFSYN_POINT_VALIDATION_H
#define RFSYN_POINT_VALIDATION_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.h"
#include "cpml.h"

static int point_validation(int margin) {
 size_t nt=margin?50:40,nr=margin?150:120,steps=margin?200:150,count=margin?5:1;
 const size_t offsets[5]={6,9,12,15,16}; int status=1;FILE *out=NULL;
 fdtd_grid_t *r=NULL,*g=NULL;cpml_t *p=NULL;double *a=NULL,*b=NULL;
 r=fdtd_grid_create(nr,nr,nr,1e-3,0.99);g=fdtd_grid_create(nt,nt,nt,1e-3,0.99);
 if(!r||!g) goto cleanup;
 p=cpml_create(g,cpml_default_params(8));a=calloc(steps*count,sizeof(double));b=calloc(steps*count,sizeof(double));
 if(!p||!a||!b) goto cleanup;
 for(size_t n=0;n<steps;n++) {
  fdtd_update_h(r);fdtd_update_h_cpml(g,p);

  double u=((double)n-40)/8,v=-2*u*exp(-u*u);
  r->ez[fdtd_index(r,nr/2,nr/2,nr/2)]+=v;g->ez[fdtd_index(g,nt/2,nt/2,nt/2)]+=v;
  fdtd_update_e(r);fdtd_update_e_cpml(g,p);
  for(size_t q=0;q<count;q++) {
   size_t off=margin?offsets[q]:12,idx=n*count+q;
   a[idx]=r->ez[fdtd_index(r,nr/2+off,nr/2,nr/2)];b[idx]=g->ez[fdtd_index(g,nt/2+off,nt/2,nt/2)];
   if(!isfinite(a[idx])||!isfinite(b[idx])) goto cleanup;
  }
 }
 out=fopen(margin?"cpml_margin_validation.csv":"cpml_validation.csv","w");if(!out) goto cleanup;
 fprintf(out,"step,t,probe_id,dist_to_pml,probe_ref,probe_test,abs_err,db_err\n");status=0;
 for(size_t q=0;q<count;q++) {
  double peak=0,worst=0,se=0,sr=0;
  for(size_t n=0;n<steps;n++) {size_t z=n*count+q;if(fabs(a[z])>peak)peak=fabs(a[z]);}
  if(!(peak>1e-12)){status=1;continue;}
  for(size_t n=0;n<steps;n++) {
   size_t z=n*count+q;double e=fabs(b[z]-a[z]);if(e>worst)worst=e;se+=e*e;sr+=a[z]*a[z];
   double dec=e>0?20*log10(e/peak):-INFINITY;
   if(fprintf(out,"%zu,%.17g,%zu,%zu,%.17g,%.17g,%.17g,%.17g\n",n,n*g->dt,q,nt-8-(nt/2+(margin?offsets[q]:12)),a[z],b[z],e,dec)<0) status=1;
  }
  printf("probe %zu: max normalized %.9g, relative L2 %.9g\n",q,worst/peak,sqrt(se/sr));

  if(worst/peak>0.01)status=1;
 }
cleanup:
 if(out&&fclose(out))status=1;free(a);free(b);cpml_destroy(p);fdtd_grid_destroy(r);fdtd_grid_destroy(g);
 printf("%s\n",status?"FAIL":"PASS");return status;
}
#endif
