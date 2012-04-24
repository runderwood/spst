# spst: spatial autocorrelation, by force

- for now, this only calculates the global moran's i for a shapefile
- don't use epsg4326 (or any other unprojected crs). use something like 3785.
- distances less than one are bumped up to one.
- did i mention this only does idw?
- and that distances less than one are bumped to one?
- so you shouldn't use an unprojected lat/lon srs.


    spst file.shp colname

- on my laptop, which is decent, it takes just under 7 seconds to build the idw and calculate global moran's i for 1000 features.
- i know that's not very fast.
- you should try doing 100000 features.
- there's also a small memory leak, but i'm not sure if it's something i'm doing wrong with geos or ogr or what. i'll get around to it.
