from environmental_raster_glc import PatchExtractor

extractor = PatchExtractor('/home/data/rasters_GLC19', size=64, verbose=True)

extractor.append('chbio_1')
extractor.append('text')
extractor.append('clc')
extractor.append('bs_top')
extractor.append('oc_top')

# extractor.add_all()

extractor.plot((43.61, 3.88))
