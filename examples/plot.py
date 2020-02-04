from environmental_raster_glc import PatchExtractor

extractor = PatchExtractor('/home/data/rasters_GLC20/archive/soilgrids/', size=256, verbose=True)

extractor.append('bdticm')
# extractor.append('text')
# extractor.append('clc')
# extractor.append('bs_top')
# extractor.append('oc_top')

# extractor.add_all()
extractor.plot((37.746517, -122.423786))
