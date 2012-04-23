#include <assert.h>
#include "sp.h"

int sp_ctxt_init(sp_ctxt* c, const char* p) {
    assert(c && c->initd != 1);
    OGRRegisterAll();
    c->initd = 1;
    c->ds = OGROpen(p, FALSE, NULL);
    if(c->ds == NULL) return -1;
    c->lnum = OGR_DS_GetLayerCount(c->ds);
    return 0;
}

int sp_ctxt_dest(sp_ctxt* c) {
    assert(c && c->initd);
    if(c->ds) OGR_DS_Destroy(c->ds);
    return 0;
}

int sp_ctxt_lstat(sp_ctxt* c, int l, sp_layerst* stat, const char* xcolnm) {
    OGRLayerH lyr = OGR_DS_GetLayer(c->ds, l);
    if(lyr == NULL) {
        return SP_ERR_NOEXIST;
    }
    int fnum = OGR_L_GetFeatureCount(lyr, TRUE);
    if(fnum <= 0) return -1;
    sp_idw* idw = calloc(sizeof(sp_idw), 1);
    assert(!sp_idw_init(idw, fnum));
    long i, j;
    int col;
    OGRFeatureH fp_i, fp_j;
    OGRGeometryH p_i, p_j;
    double X_sum = 0.0;
    double X_i, X_j;
    double* w_ij;
    for(i=0; i<fnum; i++) {
        fp_i = OGR_L_GetFeature(lyr, i);
        p_i = OGR_F_GetGeometryRef(fp_i);
        col = OGR_F_GetFieldIndex(fp_i, xcolnm);
        if(col < 0) {
            sp_idw_dest(idw);
            return SP_ERR_NOFIELD;
        }
        X_i = OGR_F_GetFieldAsDouble(fp_i, col);
        X_sum += X_i;
        for(j=0; j<fnum; j++) {
            w_ij = idw->w+i*fnum+j;
            if(i == j) {
                *w_ij = 0.0;
                continue;
            }
            fp_j = OGR_L_GetFeature(lyr, j);
            p_j = OGR_F_GetGeometryRef(fp_j);
            *w_ij = 1.0/OGR_G_Distance(p_i, p_j);
            if(isinf(*w_ij)) *w_ij = 0.0;
            OGR_F_Destroy(fp_j);
        }
        OGR_F_Destroy(fp_i);
    }

    double s0, s1, s2, I, EI, N, X_bar;
    N = (double)(fnum);
    X_bar = X_sum / (double)fnum;

    s0 = 0.0;
    s1 = 0.0;
    s2 = 0.0;

    double v_s1, v_s2, v_s2_i, v_s3, v_s4, v_s5 = 0.0;

    for(i=0; i<fnum; i++) {
        fp_i = OGR_L_GetFeature(lyr, i);
        p_i = OGR_F_GetGeometryRef(fp_i);
        col = OGR_F_GetFieldIndex(fp_i, xcolnm);
        if(col < 0) {
            sp_idw_dest(idw);
            return SP_ERR_NOFIELD;
        }
        s2 += (X_i - X_bar) * (X_i - X_bar);
        X_i = OGR_F_GetFieldAsDouble(fp_i, col);

        v_s2_i = 0.0;

        for(j=0; j<fnum; j++) {
            fp_j = OGR_L_GetFeature(lyr, j);
            p_j = OGR_F_GetGeometryRef(fp_j);
            X_j = OGR_F_GetFieldAsDouble(fp_j, col);
            w_ij = idw->w+i*fnum+j;
            w_ji = idw->w+j*fnum+i;
            s0 += *w_ij;
            s1 += (*w_ij) * (X_i - X_bar) * (X_j - X_bar);
            v_s1 += (*w_ij + *w_ji) * (*w_ij + *w_ji);
            v_s2_i += *w_ij + *w_ji;
            OGR_F_Destroy(fp_j);
        }

        v_s2 += v_s2_i * v_s2_i;

        OGR_F_Destroy(fp_i);
    }

    v_s1 *= 0.5;
    
    sp_idw_dest(idw);
    free(idw);
    idw = NULL;

    I = (N / s0) * (s1 / s2);
    EI = -1.0 / (N - 1.0);

    stat->moransi = I;
    stat->moransei = EI;

    return 0;
}

int sp_idw_init(sp_idw* idw, int n) {
    assert(!idw->w);
    idw->w = malloc(sizeof(double)*n*n);
    if(idw->w == NULL) return SP_ERR_OUTOMEM;
    idw->n = n;
    idw->sz = n*n;
    return 0;
}

int sp_idw_dest(sp_idw* idw) {
    assert(idw && idw->w);
    free(idw->w);
    idw->w = NULL;
    idw->n = 0;
    idw->sz = 0;
    return 0;
}

int sp_idw_setall(sp_idw* idw, double x) {
    assert(idw);
    int i, o, j;
    for(i=0; i<idw->n; i++) {
        for(o=0; o<idw->n; o++) {
            j = idw->n*i+o;
            *(idw->w+(idw->n*i+o)) = x;
        }
    }
    return 0;
}
