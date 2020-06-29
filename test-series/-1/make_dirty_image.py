import argparse

parser = argparse.ArgumentParser(description='Run with casa as: casa -c ' +
                                 'make_dirty_image.py <args>"')
parser.add_argument('--script', '-c', help='Name of this script', required=True)
parser.add_argument('--fitsfile', help='Path to input uvfits file', required=True)
parser.add_argument('--vis', help='Path to output measurement set (default: replace ' +
                    '.uvfits with .ms in --fitsfile)', required=False)
parser.add_argument('--imagename', help='Prename of image file (default: --fitsfile with' +
                    '.uvfits removed', required=False)
parser.add_argument('--imsize', type=int, default=2048, help='Image size (default: 2048)',
                    required=False)
parser.add_argument('--cell', default='1arcmin', help='Size of cell (default: 1arcmin)',
                    required=False)
parser.add_argument('--weighting', default='natural', help='Weighting scheme ' +
                    '(default: natural)', required=False)
parser.add_argument('--robust', default=0.5, help='Robust parameter for Briggs ' +
                    ' weighting (default: 0.5)', required=False)
args = parser.parse_args()

# Handle default paths
vis = args.fitsfile.replace('.uvfits', '.ms')
imagename = args.fitsfile.replace('.uvfits', '')
if args.vis:
    vis = args.vis
if args.imagename:
    imagename = args.imagename

# Convert uvfits to ms
print('Writing {}'.format(vis))
importuvfits(fitsfile=args.fitsfile, vis=vis)

# CLEAN
tclean(vis=vis, imagename=imagename, imsize=args.imsize, cell=args.cell,
       weighting=args.weighting, robust=args.robust, niter=0)

# Write image to fits file
imagename = imagename + '.image'
fitsimage = imagename + '.fits'
print('Writing image to {}'.format(fitsimage))
exportfits(imagename=imagename, fitsimage=fitsimage)
