#include <stdlib.h>
#include <stdio.h>
#include "sp.h"

int usage() {
    printf("Usage:\n\tspst [-w rook|ones|dist*] srcfile.shp\n");
    return 0;
}

int parse_opts(int argc, const char** argv) {
    int optflags = 0;
    return optflags;
}

int main(int argc, const char** argv) {
    if(argc <=2) {
        usage();
        goto fail;
    }
    //int optflags = parse_opts(argc, argv);
    const char* xcolnm = argv[argc-1]; 
    const char* infile_path = argv[argc-2];
    sp_ctxt* ctxt = malloc(sizeof(sp_ctxt));
    if(sp_ctxt_init(ctxt, infile_path)) {
        fprintf(stderr, "could not open source: %s\n", infile_path);
        goto fail;
    }
    fprintf(stderr, "using infile: %s (%d layer(s))\n", infile_path, ctxt->lnum);

    sp_layerst* st = malloc(sizeof(sp_layerst));
    int stres;
    if((stres = sp_ctxt_lstat(ctxt, 0, st, xcolnm))) {
        fprintf(stderr, "could not calculate stats: %d\n", stres);
        goto fail;
    }

    printf("Spatial Autocorrelation\n==========================\n");
    printf("I:\t\t%f\n", st->moransi);
    printf("E(I):\t\t%f\n", st->moransei);

    free(st);
    st = NULL;

    if(sp_ctxt_dest(ctxt)) {
        fprintf(stderr, "bad shutdown");
        goto fail;
    }
    free(ctxt);
    ctxt = NULL;

    return EXIT_SUCCESS;
    fail:
        return EXIT_FAILURE;
}

