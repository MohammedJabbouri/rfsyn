#include "material_catalog.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { color_rgb_t color; char name[MATERIAL_CATALOG_NAME_SIZE]; } color_entry_t;
struct material_catalog {
    catalog_definition_t entries[MATERIAL_CATALOG_LIMIT];
    size_t count;
    color_entry_t colors[MATERIAL_CATALOG_COLOR_LIMIT];
    size_t color_count;
};
typedef struct { const char *name, *target; } alias_t;
static const alias_t aliases[] = {
    {"vacuum","air"}, {"drywall","plasterboard"}, {"dirt","soil_medium_dry"},
    {"grass_dry_proxy","soil_very_dry"}, {"grass_moist_proxy","soil_medium_dry"},
    {"human_phantom","tissue_muscle"}, {"wood","wood_framing_spf"},
    {"brick","brick_clay"}, {"concrete","concrete_normal"}
};
static const char *canonical(const char *name) {
    if (!name) return NULL;
    for (size_t i=0;i<sizeof(aliases)/sizeof(aliases[0]);i++)
        if (!strcmp(name,aliases[i].name)) return aliases[i].target;
    return name;
}
static int text_copy(char *dst,size_t capacity,const char *src) {
    if (!src) src="";
    size_t n=strlen(src);
    if (n>=capacity) return CATALOG_INVALID;
    memcpy(dst,src,n+1);
    return CATALOG_OK;
}
static int valid_definition(const catalog_definition_t *d) {
    if (!d || !memchr(d->name,0,sizeof(d->name)) || !memchr(d->description,0,sizeof(d->description)) || !memchr(d->source,0,sizeof(d->source))) return 0;
    if (!d->name[0] || !d->source[0]) return 0;
    for (const unsigned char *p=(const unsigned char *)d->name;*p;p++)
        if (!((*p>='a' && *p<='z') || (*p>='0' && *p<='9') || *p=='_')) return 0;
    if (d->status<CATALOG_ANALYTIC || d->status>CATALOG_PENDING) return 0;
    if (!isfinite(d->min_frequency_hz) || !isfinite(d->max_frequency_hz)) return 0;
    if (!(d->min_frequency_hz==0 && d->max_frequency_hz==0) && !(d->min_frequency_hz>0 && d->max_frequency_hz>=d->min_frequency_hz)) return 0;
    if (d->model==CATALOG_UNCONFIGURED) return d->status==CATALOG_PENDING;
    if (d->status==CATALOG_PENDING) return 0;
    if (d->model==CATALOG_CONSTANT)
        return isfinite(d->eps_r) && d->eps_r>=1 && isfinite(d->sigma_s_per_m) && d->sigma_s_per_m>=0 && isfinite(d->reference_frequency_hz) && d->reference_frequency_hz>0;
    if (d->model==CATALOG_POWER_LAW)
        return isfinite(d->a) && d->a>0 && isfinite(d->b) && isfinite(d->c) && d->c>=0 && isfinite(d->d) && d->min_frequency_hz>0;
    return 0;
}
int material_catalog_find(const material_catalog_t *cat,const char *name,catalog_material_id_t *out) {
    if (!cat || !name || !out) return CATALOG_INVALID;
    name=canonical(name);
    for (size_t i=0;i<cat->count;i++) if (!strcmp(name,cat->entries[i].name)) {
        *out=(catalog_material_id_t)i; return CATALOG_OK;
    }
    return CATALOG_NOT_FOUND;
}
int material_catalog_get(const material_catalog_t *cat,catalog_material_id_t id,catalog_definition_t *out) {
    if (!cat || !out) return CATALOG_INVALID;
    if ((size_t)id>=cat->count) return CATALOG_NOT_FOUND;
    *out=cat->entries[id]; return CATALOG_OK;
}
int material_catalog_define(material_catalog_t *cat,const catalog_definition_t *d,int override_existing,catalog_material_id_t *out) {
    if (!cat || !out || !valid_definition(d)) return CATALOG_INVALID;
    if (strcmp(canonical(d->name),d->name)) return CATALOG_RESERVED;
    catalog_material_id_t id;
    if (material_catalog_find(cat,d->name,&id)==CATALOG_OK) {
        if (!strcmp(d->name,"air")) return CATALOG_RESERVED;
        if (!override_existing) return CATALOG_EXISTS;
        cat->entries[id]=*d;
    } else {
        if (cat->count==MATERIAL_CATALOG_LIMIT) return CATALOG_FULL;
        id=(catalog_material_id_t)cat->count;
        cat->entries[cat->count++]=*d;
    }
    *out=id; return CATALOG_OK;
}
int material_catalog_define_constant(material_catalog_t *cat,const char *name,double eps,double sigma,double ref,const char *description,const char *source,int override_existing,catalog_material_id_t *out) {
    catalog_definition_t d={0};
    if (!name || text_copy(d.name,sizeof(d.name),name) || text_copy(d.description,sizeof(d.description),description) || text_copy(d.source,sizeof(d.source),source?source:"User supplied; unverified")) return CATALOG_INVALID;
    d.model=CATALOG_CONSTANT; d.status=CATALOG_USER_DEFINED;
    d.eps_r=eps; d.sigma_s_per_m=sigma; d.reference_frequency_hz=ref;
    return material_catalog_define(cat,&d,override_existing,out);
}
int material_catalog_evaluate(const material_catalog_t *cat,const char *name,double frequency,int allow_proxy,catalog_properties_t *out) {
    if (!out || !isfinite(frequency) || frequency<=0) return CATALOG_INVALID;
    catalog_material_id_t id;
    int rc=material_catalog_find(cat,name,&id);
    if (rc!=CATALOG_OK) return rc;
    const catalog_definition_t *d=&cat->entries[id];
    if (d->model==CATALOG_UNCONFIGURED) return CATALOG_NEEDS_PROPERTIES;
    if (d->status==CATALOG_PROXY && !allow_proxy) return CATALOG_PROXY_DISABLED;
    if (d->min_frequency_hz>0 && (frequency<d->min_frequency_hz || frequency>d->max_frequency_hz)) return CATALOG_FREQUENCY_RANGE;
    double eps=d->eps_r, sigma=d->sigma_s_per_m;
    if (d->model==CATALOG_POWER_LAW) {
        double f=frequency/1e9;
        eps=d->a*pow(f,d->b);
        sigma=d->c==0?0:d->c*pow(f,d->d);
    }
    if (!isfinite(eps) || eps<1 || !isfinite(sigma) || sigma<0) return CATALOG_INVALID;
    *out=(catalog_properties_t){id,d->status,frequency,eps,sigma};
    return CATALOG_OK;
}
static int add_pending(material_catalog_t *cat,const char *name,const char *description) {
    catalog_definition_t d={0}; catalog_material_id_t id;
    if (text_copy(d.name,sizeof(d.name),name) || text_copy(d.description,sizeof(d.description),description) || text_copy(d.source,sizeof(d.source),"Unconfigured: documented properties required")) return CATALOG_INVALID;
    d.model=CATALOG_UNCONFIGURED; d.status=CATALOG_PENDING;
    return material_catalog_define(cat,&d,0,&id);
}
static int add_power(material_catalog_t *cat,const char *name,const char *description,catalog_status_t status,double a,double c,double exponent,double low,double high) {
    catalog_definition_t d={0}; catalog_material_id_t id;
    if (text_copy(d.name,sizeof(d.name),name) || text_copy(d.description,sizeof(d.description),description) || text_copy(d.source,sizeof(d.source),"ITU-R P.2040-3 Table 3; generic class, MathWorks reproduction")) return CATALOG_INVALID;
    d.model=CATALOG_POWER_LAW; d.status=status;
    d.a=a; d.c=c; d.d=exponent;
    d.min_frequency_hz=low; d.max_frequency_hz=high;
    return material_catalog_define(cat,&d,0,&id);
}
static int add_constant(material_catalog_t *cat,const char *name,const char *description,catalog_status_t status,const char *source,double eps,double sigma,double low,double high) {
    catalog_definition_t d={0}; catalog_material_id_t id;
    if (text_copy(d.name,sizeof(d.name),name) || text_copy(d.description,sizeof(d.description),description) || text_copy(d.source,sizeof(d.source),source)) return CATALOG_INVALID;
    d.model=CATALOG_CONSTANT; d.status=status;
    d.eps_r=eps; d.sigma_s_per_m=sigma; d.reference_frequency_hz=1e9;
    d.min_frequency_hz=low; d.max_frequency_hz=high;
    return material_catalog_define(cat,&d,0,&id);
}
material_catalog_t *material_catalog_create_default(void) {
    material_catalog_t *cat=calloc(1,sizeof(*cat));
    if (!cat) return NULL;
    if (add_constant(cat,"air","Air represented as vacuum; nonmagnetic",CATALOG_ANALYTIC,"Analytic vacuum approximation",1,0,0,0)) goto fail;
    if (add_power(cat,"wood_framing_spf","Generic wood proxy; not SPF-specific, isotropic",CATALOG_PROXY,1.99,.0047,1.0718,1e6,100e9)) goto fail;
    if (add_power(cat,"wood_pine","Generic wood proxy; not pine-specific or a living-tree model",CATALOG_PROXY,1.99,.0047,1.0718,1e6,100e9)) goto fail;
    if (add_power(cat,"wood_oak","Generic wood proxy; not oak-specific",CATALOG_PROXY,1.99,.0047,1.0718,1e6,100e9)) goto fail;
    if (add_power(cat,"brick_clay","Generic brick proxy; not clay-product-specific",CATALOG_PROXY,3.91,.0238,.16,1e9,10e9)) goto fail;
    if (add_power(cat,"brick_concrete","Generic concrete proxy; not masonry-brick-specific",CATALOG_PROXY,5.24,.0462,.7822,1e9,100e9)) goto fail;
    if (add_power(cat,"concrete_normal","Generic concrete; moisture and aggregate unspecified",CATALOG_SOURCED,5.24,.0462,.7822,1e9,100e9)) goto fail;
    if (add_power(cat,"concrete_lightweight","PROVISIONAL ordinary-concrete proxy; NOT a lightweight-concrete measurement",CATALOG_PROXY,5.24,.0462,.7822,1e9,100e9)) goto fail;
    if (add_constant(cat,"asphalt","Dry asphalt provisional GPR proxy; selected values, not validated GHz loss; 1 GHz only",CATALOG_PROXY,"GPRRental comparison table: dry asphalt eps=3, sigma=.001-.01 mS/m; selected midpoint",3,5.5e-6,1e9,1e9)) goto fail;
    if (add_pending(cat,"ceramic_tile","Ceramic tile: deferred")) goto fail;
    if (add_pending(cat,"vinyl_flooring","Vinyl flooring: deferred")) goto fail;
    if (add_pending(cat,"carpet","Carpet: deferred")) goto fail;
    if (add_pending(cat,"gravel","Gravel: deferred")) goto fail;
    if (add_constant(cat,"sand_dry","Dry sand provisional GPR proxy; selected eps=4, sigma=.01 mS/m; 1 GHz only",CATALOG_PROXY,"GPRRental comparison table: dry sand eps=3-5, sigma=.01 mS/m; not validated GHz data",4,1e-5,1e9,1e9)) goto fail;
    if (add_constant(cat,"soil_very_dry","Legacy silty-loam point approximation; mv=.07, T=23 C; 1 GHz only, pending independent verification",CATALOG_PROXY,"Prior project computation attributed to P.527-6; implementation not independently verified",4.280,.02665,1e9,1e9)) goto fail;
    if (add_constant(cat,"soil_medium_dry","Legacy silty-loam point approximation; mv=.20, T=23 C; 1 GHz only, pending independent verification",CATALOG_PROXY,"Prior project computation attributed to P.527-6; implementation not independently verified",9.896,.06674,1e9,1e9)) goto fail;
    if (add_constant(cat,"soil_wet","Legacy silty-loam point approximation; mv=.50, T=23 C; 1 GHz only, pending independent verification",CATALOG_PROXY,"Prior project computation attributed to P.527-6; implementation not independently verified",30.290,.17151,1e9,1e9)) goto fail;
    if (add_power(cat,"plasterboard","Generic plasterboard class",CATALOG_SOURCED,2.73,.0085,.9395,1e9,100e9)) goto fail;
    if (add_power(cat,"glass","Generic glass class; lower-frequency branch",CATALOG_SOURCED,6.31,.0036,1.3394,.1e9,100e9)) goto fail;
    if (add_power(cat,"ceiling_board","Generic ceiling board; lower-frequency branch",CATALOG_SOURCED,1.48,.0011,1.0750,1e9,100e9)) goto fail;
    if (add_constant(cat,"tissue_muscle","Muscle; Gabriel four-Cole-Cole fit evaluated at 1 GHz; point model only",CATALOG_SOURCED,"Gabriel 1996 AL/OE-TR-1996-0037 Appendix C; computed at 1 GHz",54.8110707802934,0.978201234619143,1e9,1e9)) goto fail;
    if (add_constant(cat,"tissue_fat","Fat (not breast fat); Gabriel four-Cole-Cole fit evaluated at 1 GHz; point model only",CATALOG_SOURCED,"Gabriel 1996 AL/OE-TR-1996-0037 Appendix C; computed at 1 GHz",5.44703776705381,0.0535028309253921,1e9,1e9)) goto fail;
    if (add_constant(cat,"tissue_skin","Dry skin (not wet skin); Gabriel four-Cole-Cole fit evaluated at 1 GHz; point model only",CATALOG_SOURCED,"Gabriel 1996 AL/OE-TR-1996-0037 Appendix C; computed at 1 GHz",40.9361354522536,0.899791175272886,1e9,1e9)) goto fail;
    return cat;
fail:
    material_catalog_destroy(cat); return NULL;
}
void material_catalog_destroy(material_catalog_t *cat) { free(cat); }
size_t material_catalog_count(const material_catalog_t *cat) { return cat?cat->count:0; }
static int same_color(color_rgb_t a,color_rgb_t b) { return a.r==b.r && a.g==b.g && a.b==b.b; }
int material_catalog_map_color(material_catalog_t *cat,color_rgb_t color,const char *name) {
    catalog_material_id_t id;
    int rc=material_catalog_find(cat,name,&id);
    if (rc!=CATALOG_OK) return rc;
    size_t slot;
    for (slot=0;slot<cat->color_count;slot++) if (same_color(cat->colors[slot].color,color)) break;
    if (slot==cat->color_count) {
        if (cat->color_count==MATERIAL_CATALOG_COLOR_LIMIT) return CATALOG_FULL;
        cat->color_count++;
    }
    cat->colors[slot].color=color;
    return text_copy(cat->colors[slot].name,sizeof(cat->colors[slot].name),cat->entries[id].name);
}
int material_catalog_lookup_color(const material_catalog_t *cat,color_rgb_t color,char *out,size_t capacity) {
    if (!cat || !out || !capacity) return CATALOG_INVALID;
    for (size_t i=0;i<cat->color_count;i++) if (same_color(cat->colors[i].color,color))
        return text_copy(out,capacity,cat->colors[i].name);
    return CATALOG_NOT_FOUND;
}
const char *material_catalog_error(int rc) {
    switch (rc) {
        case CATALOG_OK: return "success";
        case CATALOG_INVALID: return "invalid definition, frequency, or argument";
        case CATALOG_NOT_FOUND: return "unknown material";
        case CATALOG_EXISTS: return "material exists; explicit override required";
        case CATALOG_FULL: return "catalog or palette capacity reached";
        case CATALOG_NEEDS_PROPERTIES: return "material has no configured physical properties";
        case CATALOG_PROXY_DISABLED: return "generic proxy requires explicit permission";
        case CATALOG_FREQUENCY_RANGE: return "frequency outside configured model range";
        case CATALOG_RESERVED: return "reserved air or alias name; use canonical material name";
        default: return "unknown catalog error";
    }
}
