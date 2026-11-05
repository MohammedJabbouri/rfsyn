#include "../engine/fdtd/extended_test.h"

static int errors = 0;
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x); errors++; } } while (0)

static void expect_epsilon(xt_sim *s, field_component_t component,
    size_t i, size_t j, size_t k, double expected)
{
    material_coeff_t c;
    int status = materials_get_coefficients(s->m,s->g,component,i,j,k,&c);
    CHECK(status == 0);
    if (status != 0) return;
    CHECK(isfinite(c.Cb) && c.Cb > 0.0);
    CHECK(fabs(c.Ca-1.0) < 1e-14);
    double measured = s->g->dt/(s->g->dx*c.Cb*EPS0);
    if (fabs(measured-expected) > 1e-11*fmax(1.0,expected)) {
        fprintf(stderr,"coefficient mismatch component=%d (%zu,%zu,%zu): %.17g expected %.17g\n",
                (int)component,i,j,k,measured,expected);
        errors++;
    }
}

static void interface_map(unsigned axes)
{
    xt_sim s = {0};
    int status=xt_create(&s,8,6,4,1e-3,axes,0);
    CHECK(status==0); if(status!=0)return;
    status=xt_interface(&s,4,1.0,25.0);
    CHECK(status==0);
    if (status==0) {
        for(size_t i=0;i<8;i++) for(size_t j=0;j<6;j++) for(size_t k=0;k<4;k++) {
            double base=i<4?1.0:25.0;
            double tangent=(i==4 || ((axes&1u)&&i==0))?13.0:base;
            expect_epsilon(&s,FIELD_EX,i,j,k,base);
            expect_epsilon(&s,FIELD_EY,i,j,k,tangent);
            expect_epsilon(&s,FIELD_EZ,i,j,k,tangent);
        }
        if (!(axes&1u)) {
            expect_epsilon(&s,FIELD_EY,8,0,0,25.0);
            expect_epsilon(&s,FIELD_EZ,8,0,0,25.0);
        }
    }
    xt_destroy(&s);
}

static void assignment_boxes(void)
{
    xt_sim s={0};
    int status=xt_create(&s,8,6,4,1e-3,7u,0);
    CHECK(status==0);if(status!=0)return;
    material_id_t dielectric;
    CHECK(materials_register(s.m,25.0,0.0,&dielectric)==0);
    CHECK(materials_set_component(s.m,s.g,FIELD_EY,0,0,0,7.0,0.0)==0);
    CHECK(materials_set_component(s.m,s.g,FIELD_EY,4,3,2,9.0,0.0)==0);
    CHECK(materials_assign_box(s.m,s.g,4,9,3,7,2,5,dielectric)==0);
    for(size_t i=0;i<8;i++) for(size_t j=0;j<6;j++) for(size_t k=0;k<4;k++) {
        double base=(i>=4&&j>=3&&k>=2)?25.0:1.0;
        expect_epsilon(&s,FIELD_EX,i,j,k,base);
        expect_epsilon(&s,FIELD_EY,i,j,k,(i==0&&j==0&&k==0)?7.0:base);
        expect_epsilon(&s,FIELD_EZ,i,j,k,base);
    }
    CHECK(materials_assign_box(s.m,s.g,8,9,0,1,0,1,dielectric)!=0);
    expect_epsilon(&s,FIELD_EY,0,0,0,7.0);
    CHECK(materials_assign_box(s.m,s.g,0,9,0,7,0,5,dielectric)==0);
    for(size_t i=0;i<8;i++) for(size_t j=0;j<6;j++) for(size_t k=0;k<4;k++)
        for(int c=0;c<3;c++) expect_epsilon(&s,(field_component_t)c,i,j,k,25.0);
    xt_destroy(&s);
}

int main(void)
{
    interface_map(7u);
    interface_map(6u);
    interface_map(0u);
    assignment_boxes();
    if(errors)fprintf(stderr,"%d periodic material checks failed\n",errors);
    else printf("PASS: periodic box assignment, overrides, interface and seam coefficients; nonperiodic control\n");
    return errors?1:0;
}