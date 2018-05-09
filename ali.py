import alipy
import glob
import sys


sci_im = sys.argv[1]
ref_im = sys.argv[2]

identifications = alipy.ident.run(ref_im, [sci_im], visu=False)

outputshape = alipy.align.shape(ref_im)
# This is simply a tuple (width, height)... you could specify any other shape.

for id in identifications:
        if id.ok == True:

                # Variant 1, using only scipy and the simple affine transorm :
                alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True)

                # Variant 2, using geomap/gregister, correcting also for distortions :
                #alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape, makepng=False)