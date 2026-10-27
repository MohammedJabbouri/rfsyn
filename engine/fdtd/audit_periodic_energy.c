#include "extended_test.h"

static int run_audit(const char *name, double eps_r, int interface, size_t steps) {
    xt_sim s = {0};
    FILE *f = NULL;
    double *eps = NULL, *hprev = NULL;
    int rc = 1;
    const size_t N = 8, count = N*N*N, components = 3*count;
    double reference = 0.0, invariant_reference = 0.0;
    double min_ratio = INFINITY, max_ratio = -INFINITY, max_invariant_drift = 0.0;
    if (xt_create(&s,N,N,N,1e-3,7u,0) != 0) goto done;
    if ((interface ? xt_interface(&s,N/2,1.0,eps_r) : xt_medium(&s,eps_r)) != 0) goto done;
    eps = malloc(components*sizeof(*eps));
    hprev = malloc(components*sizeof(*hprev));
    if (!eps || !hprev) goto done;
    double *E[3] = {s.g->ex,s.g->ey,s.g->ez};
    double *H[3] = {s.g->hx,s.g->hy,s.g->hz};
    char path[128];
    snprintf(path,sizeof path,"periodic_audit_%s_coefficients.csv",name);
    f = fopen(path,"w");
    if (!f) goto done;
    if (fprintf(f,"i,j,k,component,epsilon_f_per_m,epsilon_relative,Ca,Cb\n") < 0) goto done;
    uint32_t seed = 12345u;
    for (size_t i=0;i<N;i++) for (size_t j=0;j<N;j++) for (size_t k=0;k<N;k++) {
        size_t idx=fdtd_index(s.g,i,j,k), q=((i*N+j)*N+k)*3;
        for (size_t a=0;a<3;a++) {
            material_coeff_t coeff;
            if (materials_get_coefficients(s.m,s.g,(field_component_t)a,i,j,k,&coeff) != 0) goto done;
            if (!isfinite(coeff.Ca) || fabs(coeff.Ca-1.0)>1e-12 || !(coeff.Cb>0.0) || !isfinite(coeff.Cb)) goto done;
            eps[q+a] = s.g->dt/(s.g->dx*coeff.Cb);
            if (!(eps[q+a]>0.0) || !isfinite(eps[q+a])) goto done;
            if (fprintf(f,"%zu,%zu,%zu,%zu,%.17g,%.17g,%.17g,%.17g\n",i,j,k,a,eps[q+a],eps[q+a]/EPS0,coeff.Ca,coeff.Cb)<0) goto done;
            seed=1664525u*seed+1013904223u;
            E[a][idx]=((double)(seed>>8)/16777216.0-0.5)*0.1;
        }
    }
    if (fclose(f)!=0) {f=NULL;goto done;} f=NULL;
    snprintf(path,sizeof path,"periodic_audit_%s_history.csv",name);
    f=fopen(path,"w");
    if (!f) goto done;
    if (fprintf(f,"completed_steps,t_e_s,t_h_minus_s,t_h_plus_s,electric_energy_j,magnetic_minus_energy_j,magnetic_plus_energy_j,raw_minus_energy_j,raw_plus_energy_j,synchronized_energy_j,modified_energy_j,raw_minus_over_initial,modified_relative_drift,max_abs_e,max_abs_h_minus,max_abs_h_plus,source_active\n")<0) goto done;
    const double volume=s.g->dx*s.g->dx*s.g->dx;
    for (size_t n=0;n<=steps;n++) {
        for (size_t i=0;i<N;i++) for (size_t j=0;j<N;j++) for (size_t k=0;k<N;k++) {
            size_t idx=fdtd_index(s.g,i,j,k),q=((i*N+j)*N+k)*3;
            for (size_t a=0;a<3;a++) hprev[q+a]=H[a][idx];
        }
        fdtd_update_h(s.g);
        double ue=0.0,um=0.0,up=0.0,us=0.0,uc=0.0,em=0.0,hm=0.0,hp=0.0;
        for (size_t i=0;i<N;i++) for (size_t j=0;j<N;j++) for (size_t k=0;k<N;k++) {
            size_t idx=fdtd_index(s.g,i,j,k),q=((i*N+j)*N+k)*3;
            for (size_t a=0;a<3;a++) {
                double e=E[a][idx], old=hprev[q+a], next=H[a][idx], avg=0.5*(old+next);
                if (!isfinite(e)||!isfinite(old)||!isfinite(next)) goto done;
                ue+=0.5*eps[q+a]*e*e;
                um+=0.5*MU0*old*old; up+=0.5*MU0*next*next;
                us+=0.5*MU0*avg*avg; uc+=0.5*MU0*old*next;
                em=fmax(em,fabs(e));hm=fmax(hm,fabs(old));hp=fmax(hp,fabs(next));
            }
        }
        ue*=volume;um*=volume;up*=volume;us*=volume;uc*=volume;
        double raw=ue+um, modified=ue+uc;
        if (!isfinite(raw)||!isfinite(modified)||!isfinite(ue+us)) goto done;
        if (n==0) {reference=raw;invariant_reference=modified;}
        if (!(reference>0.0)||!(invariant_reference>0.0)) goto done;
        double ratio=raw/reference, drift=(modified-invariant_reference)/invariant_reference;
        min_ratio=fmin(min_ratio,ratio);max_ratio=fmax(max_ratio,ratio);
        max_invariant_drift=fmax(max_invariant_drift,fabs(drift));
        if (fprintf(f,"%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,0\n",n,n*s.g->dt,(n-0.5)*s.g->dt,(n+0.5)*s.g->dt,ue,um,up,raw,ue+up,ue+us,modified,ratio,drift,em,hm,hp)<0) goto done;
        if (n<steps) fdtd_update_e_materials(s.g,s.m);
    }
    printf("%s: raw/initial min=%.17g max=%.17g; modified maximum absolute relative drift=%.17g (diagnostic only; no accuracy gate)\n",name,min_ratio,max_ratio,max_invariant_drift);
    rc=0;
done:
    if (f && fclose(f)!=0) rc=1;
    free(eps);free(hprev);xt_destroy(&s);
    if (rc) fprintf(stderr,"Periodic audit failed: %s\n",name);
    return rc;
}
int main(int argc,char **argv) {
    size_t steps;int full;
    if (xt_parse_stability(argc,argv,&steps,&full)!=0 || full) {
        fprintf(stderr,"Usage: %s [--steps 100000] (10000..10000000); --full unsupported\n",argv[0]);return 2;
    }
    int rc=run_audit("vacuum_3d",1.0,0,steps);
    rc|=run_audit("interface25_3d",25.0,1,steps);
    return rc?1:0;
}