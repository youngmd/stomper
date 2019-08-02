#!/usr/bin/env python
import sys
import os
from astropy.io.fits import getheader
from astropy.table import Table
from astroquery.skyview import SkyView
from astropy import units as u
from astropy import coordinates
import numpy as np
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import getopt
import footprintfinder
from kapteyn import maputils

def get_matches(box):
    ngc_table = Table.read('ngc2000.fits', format='fits')
    ic_table = Table.read('ic2000.fits', format='fits')
    matches = []
    for ngc in ngc_table:
        if box[0] < ngc['RA'] < box[1] and box[2] < ngc['DEC'] < box[3]:
            # c2 = coordinates.SkyCoord(ngc['RA'], ngc['DEC'], unit=('deg', 'deg'), frame='fk5')
            # sep = c1.separation(c2)
            # if sep < radius:
            print "Match found!  NGC %s" % ngc['NGCNUM']
            matches.append(ngc)

    # for ic in ic_table:
    #     if box[0] < ic['RA'] < box[1] and box[2] < ic['DEC'] < box[3]:
    #     # c2 = coordinates.SkyCoord(ic['RA'], ic['DEC'], unit=('deg', 'deg'), frame='fk5')
    #     # sep = c1.separation(c2)
    #     # if sep < radius:
    #         print "Match found!  IC %s" % ic['ICNUM']
    #         matches.append(ic)
    return matches

def process_opts(opts):
    circles = []
    outfile = 'test.png'

    for o, a in opts:
        if o in ("-c", "--circle"):
            tmp = a.split(':')
            print tmp
            for vals in tmp:
                val = vals.split(',')
                print val
                name = val[0]
                rad = float(val[1])
                circles.append([name, rad])
        if o in ("-o"):
            outfile = a
    return [circles, outfile]


def process_args(args):
    footprints = []
    centers = []

    for infile in args:
        fname = os.path.basename(infile)

        header = getheader(infile)

        flip = False
        try:
            if header['CTYPE1'].startswith('DEC'):
                flip = True
        except:
            pass


        fpf_args = "-r " + infile
        fpf_out = footprintfinder.main(fpf_args)

        print fpf_out
        try:
            if flip:
                ra_cen = fpf_out[4]
                dec_cen = fpf_out[3]
            else:
                ra_cen = fpf_out[3]
                dec_cen = fpf_out[4]
            polygon = fpf_out[-1]
        except:
            print "Failed reading fpf output.  Exiting"
            sys.exit(1)

        coord_strings = []
        if polygon.startswith('Union'):
            polygon = polygon[polygon.find("(") + 1:polygon.find(")")]
            p = polygon.split('Polygon')
            for pg in p:
                if pg == '':
                    continue
                else:
                    coord_strings.append(pg)
        else:
            coord_strings.append(polygon.split('Polygon')[1])

        for c in coord_strings:
            words = c.split()
            footprint = []
            try:
                float(words[0])
                start = 0
            except:
                start = 1
            for i in range(start, len(words), 2):
                if flip:
                    footprint.append([float(words[i + 1]), float(words[i])])
                else:
                    footprint.append([float(words[i]), float(words[i + 1])])

            footprints.append(footprint)

        centers.append([ra_cen, dec_cen])

    return [footprints, centers]


def find_matches(centers, footprints):
    centers = np.array(centers)
    print centers[:, 0]
    ra_cen = np.mean(centers[:, 0])
    dec_cen = np.mean(centers[:, 1])

    print ra_cen, dec_cen

    c = coordinates.SkyCoord(ra_cen, dec_cen, unit=('deg', 'deg'), frame='fk5')

    max_sep = 0.0 * u.deg
    min_ra = ra_cen
    max_ra = ra_cen
    min_dec = dec_cen
    max_dec = dec_cen

    for foot in footprints:
        for f in foot:
            if f[0] < min_ra:
                min_ra = f[0]
            if f[0] > max_ra:
                max_ra = f[0]
            if f[1] < min_dec:
                min_dec = f[1]
            if f[1] > max_dec:
                max_dec = f[1]
            fc = coordinates.SkyCoord(f[0], f[1], unit=('deg', 'deg'), frame='fk5')
            sep = c.separation(fc)
            print sep
            if sep > max_sep:
                max_sep = sep

    print "using radius %s" % max_sep
    coverage_box = [min_ra, max_ra, min_dec, max_dec]
    matches = get_matches(coverage_box)

    return [c, max_sep, matches, max_ra, min_ra]


def plot_images(c, max_sep, max_ra, min_ra, footprints, matches, circles, outfile):
    paths = SkyView.get_images(position=c, survey=['DSS'], pixels=2000, radius=max_sep * 1.8)

    for path in paths:
        path[0].writeto('/tmp/tmp.fits', clobber=True)
        f = maputils.FITSimage('/tmp/tmp.fits')
        fig = plt.figure()
        frame = fig.add_subplot(1, 1, 1)

        annim = f.Annotatedimage(frame, cmap='gist_yarg')
        annim.Image()
        grat = annim.Graticule()
        grat.setp_gratline(visible=False)

        first = True
        for foot in footprints:
            if first:
                ls = 'solid'
                first = False
            else:
                ls = 'solid'
            lons = []
            lats = []
            for toe in foot:
                lons.append(toe[0])
                lats.append(toe[1])
            annim.Skypolygon(prescription=None, lons=lons, lats=lats, fill=False, alpha=0.5, ls=ls)

        for match in matches:
            try:
                s = 'NGC' + str(match['NGCNUM'])
            except:
                continue
            pos = "%s deg %s deg" % (match['RA'], match['DEC'])
            rad = match['RADIUS']
            fs = 6
            for circle in circles:
                if circle[0] == s:
                    print "placing circle!"
                    print circle[1]
                    print pos
                    annim.Skypolygon('ellipse', cpos=pos, major=2 * circle[1], minor=2 * circle[1], units='arcmin',
                                     ls='dotted', fc='m', fill=False)
                    fs = 8
            ra_size = (max_ra - min_ra) / 2000.0
            pixradsize = rad / ra_size
            print pixradsize
            x, y = annim.topixel(match['RA'], match['DEC'])
            y2 = y - 0.3 * pixradsize
            # x2 = x + pixradsize / 2.0
            plt.annotate(s, xy=(x, y), xytext=(x, y2), fontsize=fs, path_effects=[PathEffects.withStroke(linewidth=1,
                                                                                                         foreground="w")])
        as_ruler = annim.Ruler(x1=100, y1=1900, rulersize=180, lambda0=0.0,
                               step=30, units='arcsec', addangle=90, size=0)
        as_ruler.tickdy = -10
        annim.Ruler(x1=100, y1=1900, rulersize=3.1, lambda0=0.0,
                    step=1, units='arcmin', addangle=90, size=5)
        annim.Ruler(x1=100, y1=1900, rulersize=0.05, lambda0=0.0,
                                step=0.025239, units='deg', fliplabelside=True, fun=lambda x: x / 0.025239 * 10,
                                fmt="%d kpc", addangle=90, size=4)
        annim.plot()
        annim.interact_toolbarinfo()
        annim.interact_imagecolors()
        annim.interact_writepos(wcsfmt="%f", zfmt=None, pixfmt=None, hmsdms=False)

        plt.savefig(outfile)


def main():
    opts, args = getopt.getopt(sys.argv[1:], "hc:o:f", ["help", "circle=", "out="])
    circles, outfile = process_opts(opts)
    footprints, centers = process_args(args)
    center, max_sep, matches, max_ra, min_ra = find_matches(centers, footprints)

    plot_images(center, max_sep, max_ra, min_ra, footprints, matches, circles, outfile)

    print "All done!"


if __name__ == "__main__":
    main()