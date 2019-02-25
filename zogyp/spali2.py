import subprocess
from multiprocessing import Pool
import sys
import os
import inspect
import argparse
import itertools

from astropy.io import fits
from astropy.table import Table
import numpy as np
from scipy import linalg, interpolate
from scipy.ndimage import interpolation, map_coordinates
from scipy.spatial import distance


#Use the correct call for sextractor
try:
    sex = 'sextractor'
    subprocess.call(sex, stderr=subprocess.DEVNULL)
except:
   sex = 'sex'
   subprocess.call(sex, stderr=subprocess.DEVNULL) 


root_config = os.path.dirname(os.path.realpath(__file__))+'/configfls'

######################################################################################
class AffineTransform:
    """
    Represents an affine transformation consisting of rotation, isotropic
    scaling, and shift. [x', y'] = [[a -b], [b a]] * [x, y] + [c d]
    """
    def __init__(self, v):
        self.v = np.asarray(v)

    def inverse(self):
        """
        Returns the inverse transform
        """
        # To represent affine transformations with matrices,
        # we can use homogeneous coordinates.
        homo = np.array([[self.v[0], -self.v[1], self.v[2]],
                        [self.v[1],  self.v[0], self.v[3]],
                        [0.0, 0.0, 1.0]])
        inv = linalg.inv(homo)
        return AffineTransform((inv[0, 0], inv[1, 0], inv[0, 2], inv[1, 2]))

    def matrix_form(self):
        """
        Special output for scipy.ndimage.interpolation.affine_transform
        """
        return (np.array([[self.v[0], -self.v[1]],
                         [self.v[1], self.v[0]]]), self.v[2:4])

    def apply_transform(self, xy):
        """
        Applies the transform to an array of x, y points
        """
        xy = np.asarray(xy)
        # Can consistently refer to x and y as xy[0] and xy[1] if xy is
        # 2D (1D coords) or 3D (2D coords) if we transpose the 2D case of xy
        if xy.ndim == 2:
            xy = xy.T
        xn = self.v[0]*xy[0] - self.v[1]*xy[1] + self.v[2]
        yn = self.v[1]*xy[0] + self.v[0]*xy[1] + self.v[3]
        if xy.ndim == 2:
            return np.column_stack((xn, yn))
        return np.stack((xn, yn))
#########################################################################################

#########################################################################################
def find_spline_transform(a_t, sm, tm ):
    """
    Determine the residual `x` and `y` offsets between matched coordinates
    after affine transformation and fit 2D spline surfaces to describe the
    spatially-varying correction to be applied.
    """
    spline_order = 3

    # Get the source, after affine transformation, and template coordinates
    source_coo = a_t.apply_transform(get_det_coords(sm))
    template_coo = get_det_coords(tm)
    # Create splines describing the residual offsets in x and y left over
    # after the affine transformation
    kx = ky = spline_order
    sbs_x = interpolate.SmoothBivariateSpline(template_coo[:, 0],
            template_coo[:, 1],(template_coo[:, 0]- source_coo[:, 0]), 
            kx=kx, ky=ky)

    sbs_y = interpolate.SmoothBivariateSpline(template_coo[:, 0],
            template_coo[:, 1],(template_coo[:, 1]- source_coo[:, 1]),
            kx=kx, ky=ky)

    return(sbs_x, sbs_y)
#########################################################################################

#########################################################################################
def match_stars(transform, s_coo, t_coo, t_cat, s_cat, minmatchdist=5):
    source_coo_trans = transform.apply_transform(s_coo)

    dists = distance.cdist(source_coo_trans, t_coo)
    mindists = np.min(dists, axis=1)
    passed = mindists <= minmatchdist
    sorted_idx = np.argsort(dists[passed, :])

    nmatched = np.sum(passed)
    source_matchstars = s_cat[passed]
    template_matchstars = t_cat[sorted_idx[:, 0]]

    return(nmatched, source_matchstars, template_matchstars)
#########################################################################################


#########################################################################################
def find_affine_transform(t_quad, s_quad,s_coo, t_coo, t_cat, s_cat, 
                          minmatchdist=5, minnmatch=50, maxcands=10, 
                          minquaddist=0.005):

    """
    Use the quadlist hashes to determine an initial guess at an affine
    transformation and determine matched detections lists. Then refine
    the transformation using the matched detection lists.
    """
    template_hash = np.array([q[1] for q in t_quad])
    source_hash = np.array([q[1] for q in s_quad])

    dists = distance.cdist(template_hash, source_hash)
    minddist_idx = np.argmin(dists, axis=0)
    mindist = np.min(dists, axis=0)
    best = np.argsort(mindist)

    # Use best initial guess at transformation to get list of matched stars
    for i in range(min(maxcands, len(best))):
        bi = best[i]
        template_quad = t_quad[minddist_idx[bi]]
        source_quad = s_quad[bi]
        # Get a quick (exact) transformation guess
        # using first two detections
        transform = calc_affine_transform(source_quad[0][:2],
                                        template_quad[0][:2])
        dist = mindist[bi]
        print(dist)
        passed = False
        if dist < minquaddist:
            nmatch, source_matchdets, template_matchdets = \
              match_stars(transform, s_coo, t_coo, t_cat, s_cat,  minmatchdist=minmatchdist)
            if nmatch > minnmatch:
                passed = True
                break

    if passed == False:
        raise ValueError('Not enough matching stars!')


    else:
        # Refine the transformation using the matched detections
        source_match_coo = get_det_coords(source_matchdets)
        template_match_coo = get_det_coords(template_matchdets)
        transform = calc_affine_transform(source_match_coo,
                                        template_match_coo)
        # Store the final matched detection tables and transform
        nmatch, source_matchdets, template_matchdets = \
        match_stars(transform, s_coo, t_coo, t_cat, s_cat, minmatchdist=minmatchdist)
        affine_transform = transform

        return(affine_transform, nmatch, source_matchdets, template_matchdets)

#########################################################################################

#########################################################################################
def calc_affine_transform(source_coo, template_coo):
    """
    Calculates the affine transformation
    """
    Len = len(source_coo)
    template_matrix = template_coo.ravel()
    source_matrix = np.zeros((Len*2, 4))
    source_matrix[::2, :2] = np.column_stack((source_coo[:, 0],
                                           - source_coo[:, 1]))
    source_matrix[1::2, :2] = np.column_stack((source_coo[:, 1],
                                             source_coo[:, 0]))
    source_matrix[:, 2] = np.tile([1, 0], Len)
    source_matrix[:, 3] = np.tile([0, 1], Len)

    if Len == 2:
        transform = linalg.solve(source_matrix, template_matrix)
    else:
        transform = linalg.lstsq(source_matrix, template_matrix)[0]

    return(AffineTransform(transform))
#########################################################################################

#########################################################################################
def run_sex(img_name):
    """
    Runs Sextractor
    """
    catname = (img_name+'.cat')
    print('Finding soure locations of '+img_name)
    subprocess.call([sex, '-c', root_config+'/spali.sex', img_name, '-CATALOG_NAME', catname,
                     '-PARAMETERS_NAME', root_config+'/spali.param', 
                     '-STARNNW_NAME', root_config+'/default.nnw',
                     '-FILTER_NAME', root_config+'/default.conv', 
                      ], stderr=subprocess.DEVNULL)

    return (catname)
#########################################################################################

#########################################################################################
def get_ndets(cat1, cat2):
    """
    Determines number of tections
    """
    Acat = Table.read(cat1, format='ascii.sextractor')
    Bcat = Table.read(cat2, format='ascii.sextractor')
    Sml_count = min(len(Acat), len(Bcat))
    ndets = int(0.5*Sml_count)
    return(ndets)
 
#########################################################################################


#########################################################################################
def trim_cat(ndets, cat, minfwhm=2, maxflag=4):
    """
    Trim a detection catalogue based on some SExtractor values.
    Sort this by the brightest objects then cut to the to
    """
    X = 'X_IMAGE'
    Y = 'Y_IMAGE'
    FLUX = 'FLUX_BEST'
    FWHM = 'FWHM_IMAGE'
    FLAG = 'FLAGS'
    COLUMNS = [X, Y, FLUX, FWHM, FLAG]

    redcat = Table.read(cat, format='ascii.sextractor')
    redcat = redcat[COLUMNS]
    redcat = redcat[(redcat[FWHM] > minfwhm)
             & (redcat[FLAG] < maxflag)]
    redcat.sort(FLUX)
    redcat.reverse()
    redcat = redcat[:ndets]
    return(redcat)
#########################################################################################


#########################################################################################
def get_det_coords(cat):
    return cat['Y_IMAGE', 'X_IMAGE'].as_array().view((float, 2))
#########################################################################################

#########################################################################################

def quad(combo, dist):
    """
    Create a hash from a combination of four stars (a "quad").
    References
     ----------
    Based on the algorithm of [L10]_.
    [L10] Lang, D. et al. "Astrometry.net: Blind astrometric
    calibration of arbitrary astronomical images", AJ, 2010.
    """
    max_dist_idx = np.argmax(dist)
    orders = [(0, 1, 2, 3),
              (0, 2, 1, 3),
              (0, 3, 1, 2),
              (1, 2, 0, 3),
              (1, 3, 0, 2),
              (2, 3, 0, 1)]
    order = orders[max_dist_idx]
    combo = combo[order, :]
    # Look for matrix transform [[a -b], [b a]] + [c d]
    # that brings A and B to 00 11 :
    x = combo[1, 0] - combo[0, 0]
    y = combo[1, 1] - combo[0, 1]
    b = (x-y) / (x**2 + y**2)
    a = (1/x) * (1 + b*y)
    c = b*combo[0, 1] - a*combo[0, 0]
    d = -(b*combo[0, 0] + a*combo[0, 1])

    t = AffineTransform((a, b, c, d))
    (xC, yC) = t.apply_transform((combo[2, 0], combo[2, 1])).ravel()
    (xD, yD) = t.apply_transform((combo[3, 0], combo[3, 1])).ravel()
    _hash = (xC, yC, xD, yD)
    # Break symmetries if needed
    testa = xC > xD
    testb = xC + xD > 1
    if testa:
        if testb:
            _hash = (1.0-xC, 1.0-yC, 1.0-xD, 1.0-yD)
            order = (1, 0, 2, 3)
        else:
            _hash = (xD, yD, xC, yC)
            order = (0, 1, 3, 2)
    elif testb:
        _hash = (1.0-xD, 1.0-yD, 1.0-xC, 1.0-yC)
        order = (1, 0, 3, 2)
    else:
        order = (0, 1, 2, 3)

    return combo[order, :], _hash

#########################################################################################

#########################################################################################

def make_quadlist(cat, nquaddets=20, minquadsep=50):
    """
    Create a list of hashes for "quads" of the brightest sources
    in the detection catalogue.
    """
    coo = get_det_coords(cat)

    quadlist = []
    quad_idxs = itertools.combinations(range(nquaddets), 4)
    for quad_idx in quad_idxs:
        combo = coo[quad_idx, :]
        dists = distance.pdist(combo)
        if np.min(dists) > minquadsep:
            quadlist.append(quad(combo, dists))
    return(quadlist, coo)
#########################################################################################


#########################################################################################
def spline_transform(sbs_x, sbs_y, xy, relative=False):
    """
    Returns the relative shift of xy coordinates if relative is True,
    otherwise return the value of the transformed coordinates
    """
    x0 = xy[0]
    y0 = xy[1]
    if relative is True:
        x0 = y0 = 0

    return(np.array((x0 - sbs_x.ev(xy[0], xy[1]), y0 - sbs_y.ev(xy[0], xy[1]))))
#########################################################################################


#########################################################################################
def align(sci, a_t, sbs_x, sbs_y, name, data, spline):
    """
    Perform the alignment and write the transformed source
    file.
    """
    source_fits = fits.open(sci)
    source_data = fits.getdata(sci)
    shape = source_data.shape
    if spline == True:
        def final_transform(xy):
            return ((a_t.inverse().apply_transform(xy))+spline_transform(sbs_x, sbs_y, xy, relative=True))
    else:
        def final_transform(xy):
            return (a_t.inverse().apply_transform(xy))

    yy, xx = np.meshgrid(np.arange(shape[1]), np.arange(shape[0]))
    spline_coords_shift = final_transform((xx, yy))
    source_data_transform = map_coordinates(source_data,spline_coords_shift)

    source_fits[0].data = source_data_transform
    source_fits[0].header['COMMENT'] = 'Aligned using spali2 from ZiP'
    source_fits.writeto(name, overwrite=True)
    if data == True:
        return(name, source_data_transform)
    else:
        return(name)
#########################################################################################


#########################################################################################
def spalipy(sci, ref, name = 'aligned.fits', data = False, spline = True):
    """
    Credit to Joe Lyman. This is a modded version of his algorithm to 
    operate quickly with the ZOGY subtraction
    """

    p = Pool(2)
    src_cat, ref_cat = p.map(run_sex, [(sci), (ref)])
    p.close()

    ndets = get_ndets(src_cat, ref_cat)

    p = Pool(2)
    s_cat, r_cat = p.starmap(trim_cat, [(ndets, src_cat), (ndets, ref_cat)])
    p.close()

    p = Pool(2)
    #quad_coo contains the quads and coordinate catalog, respectively
    s_quad_coo, r_quad_coo = p.map(make_quadlist, [(s_cat), (r_cat)])
    p.close()

    #at = affine_transform
    at, nmatch, source_matchdets, template_matchdets = find_affine_transform(r_quad_coo[0], s_quad_coo[0], s_quad_coo[1], 
                                                                             r_quad_coo[1], r_cat, s_cat)

    if spline ==True:
        sbs_x, sbs_y = find_spline_transform(at, source_matchdets, template_matchdets)
    else:
        sbs_x = 0
        sbs_y = 0

    subprocess.call(['rm', src_cat, ref_cat])

    if data == True:
        OUT_NAME, data_a = align(sci, at, sbs_x, sbs_y, name, data, spline)
        return(OUT_NAME, data_a)
    else:	
        OUT_NAME = align(sci, at, sbs_x, sbs_y, name, data, spline)
        return(OUT_NAME)

#########################################################################################
