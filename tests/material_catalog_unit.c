#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../engine/world/material_catalog.h"
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"FAIL at line %d: %s\n",__LINE__,#c); goto cleanup; } } while (0)
static int near(double a,double b) { return isfinite(a) && isfinite(b) && fabs(a-b)<=1e-10*fmax(1,fmax(fabs(a),fabs(b))); }
int main(void) {
    int result=1;
    material_catalog_t *cat=material_catalog_create_default();
    catalog_material_id_t id,other,old;
    catalog_properties_t p,q;
    catalog_definition_t d;
    char name[MATERIAL_CATALOG_NAME_SIZE];
    const char *names[]={"air","wood_framing_spf","wood_pine","wood_oak","brick_clay","brick_concrete","concrete_normal","concrete_lightweight","asphalt","ceramic_tile","vinyl_flooring","carpet","gravel","sand_dry","soil_very_dry","soil_medium_dry","soil_wet","plasterboard","glass","ceiling_board","tissue_muscle","tissue_fat","tissue_skin"};
    CHECK(cat!=NULL);
    CHECK(material_catalog_count(cat)==23);
    for (size_t i=0;i<23;i++) {
        CHECK(material_catalog_find(cat,names[i],&id)==CATALOG_OK && id==i);
        CHECK(material_catalog_get(cat,id,&d)==CATALOG_OK && !strcmp(d.name,names[i]));
    }
    CHECK(material_catalog_evaluate(cat,"vacuum",915e6,0,&p)==CATALOG_OK && p.id==0 && p.eps_r==1 && p.sigma_s_per_m==0);
    CHECK(material_catalog_evaluate(cat,"concrete",2.4e9,0,&p)==CATALOG_OK);
    CHECK(near(p.eps_r,5.24) && near(p.sigma_s_per_m,.0462*pow(2.4,.7822)));
    CHECK(material_catalog_evaluate(cat,"concrete",915e6,0,&p)==CATALOG_FREQUENCY_RANGE);
    CHECK(material_catalog_evaluate(cat,"wood_oak",1e9,0,&p)==CATALOG_PROXY_DISABLED);
    CHECK(material_catalog_evaluate(cat,"wood_oak",1e9,1,&p)==CATALOG_OK);
    CHECK(material_catalog_evaluate(cat,"concrete_lightweight",1e9,0,&p)==CATALOG_PROXY_DISABLED);
    CHECK(material_catalog_evaluate(cat,"concrete_lightweight",1e9,1,&p)==CATALOG_OK && p.status==CATALOG_PROXY);
    CHECK(material_catalog_evaluate(cat,"asphalt",1e9,0,&p)==CATALOG_PROXY_DISABLED);
    CHECK(material_catalog_evaluate(cat,"asphalt",1e9,1,&p)==CATALOG_OK && near(p.eps_r,3) && near(p.sigma_s_per_m,5.5e-6));
    CHECK(material_catalog_evaluate(cat,"sand_dry",1e9,1,&p)==CATALOG_OK && near(p.eps_r,4) && near(p.sigma_s_per_m,1e-5));
    CHECK(material_catalog_evaluate(cat,"asphalt",2.4e9,1,&p)==CATALOG_FREQUENCY_RANGE);
    CHECK(material_catalog_evaluate(cat,"sand_dry",915e6,1,&p)==CATALOG_FREQUENCY_RANGE);
    const char *pending[]={"ceramic_tile","vinyl_flooring","carpet","gravel"};
    for (size_t i=0;i<4;i++) CHECK(material_catalog_evaluate(cat,pending[i],1e9,1,&p)==CATALOG_NEEDS_PROPERTIES);
    const char *soils[]={"soil_very_dry","soil_medium_dry","soil_wet"};
    const double soil_eps[]={4.280,9.896,30.290};
    const double soil_sigma[]={.02665,.06674,.17151};
    for (size_t i=0;i<3;i++) {
        CHECK(material_catalog_evaluate(cat,soils[i],1e9,0,&p)==CATALOG_PROXY_DISABLED);
        CHECK(material_catalog_evaluate(cat,soils[i],1e9,1,&p)==CATALOG_OK);
        CHECK(near(p.eps_r,soil_eps[i]) && near(p.sigma_s_per_m,soil_sigma[i]));
        CHECK(material_catalog_evaluate(cat,soils[i],2.4e9,1,&p)==CATALOG_FREQUENCY_RANGE);
    }
    CHECK(material_catalog_find(cat,"dirt",&id)==CATALOG_OK);
    CHECK(material_catalog_find(cat,"soil_medium_dry",&other)==CATALOG_OK && id==other);
    CHECK(material_catalog_evaluate(cat,"tissue_muscle",1e9,0,&p)==CATALOG_OK && p.status==CATALOG_SOURCED);
    CHECK(near(p.eps_r,54.8110707802934) && near(p.sigma_s_per_m,0.978201234619143));
    CHECK(material_catalog_evaluate(cat,"tissue_muscle",2.4e9,0,&p)==CATALOG_FREQUENCY_RANGE);
    CHECK(material_catalog_evaluate(cat,"tissue_fat",1e9,0,&p)==CATALOG_OK && p.status==CATALOG_SOURCED);
    CHECK(near(p.eps_r,5.44703776705381) && near(p.sigma_s_per_m,0.0535028309253921));
    CHECK(material_catalog_evaluate(cat,"tissue_fat",2.4e9,0,&p)==CATALOG_FREQUENCY_RANGE);
    CHECK(material_catalog_evaluate(cat,"tissue_skin",1e9,0,&p)==CATALOG_OK && p.status==CATALOG_SOURCED);
    CHECK(near(p.eps_r,40.9361354522536) && near(p.sigma_s_per_m,0.899791175272886));
    CHECK(material_catalog_evaluate(cat,"tissue_skin",2.4e9,0,&p)==CATALOG_FREQUENCY_RANGE);
    CHECK(material_catalog_evaluate(cat,"human_phantom",1e9,0,&p)==CATALOG_OK);
    CHECK(material_catalog_evaluate(cat,"tissue_muscle",1e9,0,&q)==CATALOG_OK && p.id==q.id);
    CHECK(material_catalog_map_color(cat,(color_rgb_t){128,128,128},"concrete")==CATALOG_OK);
    CHECK(material_catalog_lookup_color(cat,(color_rgb_t){128,128,128},name,sizeof(name))==CATALOG_OK && !strcmp(name,"concrete_normal"));
    CHECK(material_catalog_map_color(cat,(color_rgb_t){128,128,128},"air")==CATALOG_OK);
    CHECK(material_catalog_lookup_color(cat,(color_rgb_t){128,128,128},name,sizeof(name))==CATALOG_OK && !strcmp(name,"air"));
    CHECK(material_catalog_map_color(cat,(color_rgb_t){1,2,3},"unknown")==CATALOG_NOT_FOUND);
    CHECK(material_catalog_lookup_color(cat,(color_rgb_t){1,2,3},name,sizeof(name))==CATALOG_NOT_FOUND);
    CHECK(material_catalog_define_constant(cat,"my_wall",4,.01,1e9,"Test","Test",0,&id)==CATALOG_OK);
    CHECK(material_catalog_define_constant(cat,"my_wall",5,.02,1e9,"Test","Test",0,&other)==CATALOG_EXISTS);
    CHECK(material_catalog_define_constant(cat,"my_wall",5,.02,1e9,"Test","Test",1,&other)==CATALOG_OK && other==id);
    CHECK(material_catalog_find(cat,"asphalt",&old)==CATALOG_OK);
    CHECK(material_catalog_define_constant(cat,"asphalt",4,.02,1e9,"Test only","Test",1,&id)==CATALOG_OK && id==old);
    CHECK(material_catalog_evaluate(cat,"asphalt",1e9,0,&p)==CATALOG_OK && p.status==CATALOG_USER_DEFINED);
    CHECK(material_catalog_define_constant(cat,"air",2,0,1e9,"Test","Test",1,&id)==CATALOG_RESERVED);
    CHECK(material_catalog_define_constant(cat,"dirt",2,0,1e9,"Test","Test",1,&id)==CATALOG_RESERVED);
    CHECK(material_catalog_define_constant(cat,"bad",NAN,0,1e9,"Test","Test",0,&id)==CATALOG_INVALID);
    CHECK(material_catalog_define_constant(cat,"bad",2,-1,1e9,"Test","Test",0,&id)==CATALOG_INVALID);
    CHECK(material_catalog_define_constant(cat,"bad",2,0,1e9,"Test","Test",0,NULL)==CATALOG_INVALID);
    CHECK(material_catalog_evaluate(cat,"air",INFINITY,0,&p)==CATALOG_INVALID);
    CHECK(material_catalog_evaluate(cat,"air",1e9,0,NULL)==CATALOG_INVALID);
    CHECK(material_catalog_find(NULL,"air",&id)==CATALOG_INVALID);
    while (material_catalog_count(cat)<MATERIAL_CATALOG_LIMIT) {
        snprintf(name,sizeof(name),"test_%zu",material_catalog_count(cat));
        CHECK(material_catalog_define_constant(cat,name,2,0,1e9,"Test","Test",0,&id)==CATALOG_OK);
    }
    CHECK(id==255);
    CHECK(material_catalog_define_constant(cat,"overflow",2,0,1e9,"Test","Test",0,&id)==CATALOG_FULL);
    puts("material_catalog_unit: PASS"); result=0;
cleanup:
    material_catalog_destroy(cat); return result;
}