#ifndef SP_H
#define SP_H 1

#include <ogr_api.h>

#define SP_ERR_BADGEOM 1
#define SP_ERR_OUTOMEM 2
#define SP_ERR_NOEXIST 4
#define SP_ERR_NOFIELD 8

typedef struct sp_ctxt_s {
    int initd;
    int lnum;
    OGRDataSourceH* ds;
} sp_ctxt;

typedef struct sp_idw_s {
    double* w;
    int n;
    size_t sz;
} sp_idw;

typedef struct sp_layerst_s {
    double moransi;
    double moransei;
} sp_layerst;

double sp_dist(double x0, double y0, double x1, double y1);

int sp_ctxt_init(sp_ctxt* c, const char* p);
int sp_ctxt_dest(sp_ctxt* c);
int sp_ctxt_lstat(sp_ctxt* c, int l, sp_layerst* s, const char* xcolnm);

int sp_idw_init(sp_idw* idw, int n);
int sp_idw_dest(sp_idw* idw);
int sp_idw_setall(sp_idw* idw, double x);

#endif
