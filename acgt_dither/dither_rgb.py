# Original by Peter Cock from https://raw.github.com/peterjc/picobio/master/acgt_dither/dither_rgb.py
# Modified by Aaron Hardin

import os
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from PIL import Image

from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Drawing, String, Line, Rect, Wedge
from reportlab.lib.units import mm, cm, inch
from reportlab.lib import colors

#These are hueristic, and currently slighlty squash the image vertically,
#so h_scale should be a little smaller or the v_scale a little bigger.
h_scale = 0.140 * cm #per bp
v_scale = 0.125 * cm #per bp

def run(im, seq, pdf_file, main_caption):
    data = np.array(im) #.reshape(shape_rgb, order="F") 
    shape = data.shape[0:2]

    pixels = np.product(shape)
    print "Have %i base pairs, shape %r, and %i pixels" % (len(seq), shape, pixels)
    
    assert pixels <= len(seq)
    assert 0 <= data.min() <= data.max() <= 255

    #Open PDF
    width, height = page_size = h_scale * shape[1], v_scale * shape[0]
    print "Creating %s, %i by %i mm" % (pdf_file, width / mm, height / mm)
    c = canvas.Canvas(pdf_file, page_size)
    c.setTitle(main_caption)
    d = Drawing(*page_size)
    base = 0
    r_max = float(data[:,:,0].max())
    g_max = float(data[:,:,1].max())
    b_max = float(data[:,:,2].max())
    #r_mean = float(data[:,:,0].mean())
    #g_mean = float(data[:,:,1].mean())
    #b_mean = float(data[:,:,2].mean())
    #r = Rect(0,0,width,height,fillColor=colors.Color(r_mean,g_mean,b_mean))
    #d.add(r)
    for row in range(shape[0]):
        for col in range(shape[1]):
            #Ignore any alpha channel
            r, g, b = data[row, col, 0:3]
            color = colors.Color(r / r_max, g / g_max, b / b_max)
            s = String((col + 0.5) * h_scale, (shape[0]-row) * v_scale,
                       seq[base], fillColor = color,
                       fontSize = 4, textAnchor = "middle")
            d.add(s)
            base += 1
    renderPDF.draw(d, c, 0, 0)
    c.showPage()
    c.save()

#print "A0: Suggest width %i, height %i pixels" % (841 * mm / h_scale, 1189 * mm / v_scale)
#--> A0: Suggest width 600, height 951 pixels
#print "A1: Suggest width %i, height %i pixels" % (594 * mm / h_scale, 841 * mm / v_scale)
#--> Suggest width 424, height 672 pixels
#print "A2: Suggest width %i, height %i pixels" % (420 * mm / h_scale, 594 * mm / v_scale)
#--> A2: Suggest width 300, height 475 pixels
#print "A3: Suggest width %i, height %i pixels" % (297 * mm / h_scale, 420 * mm / v_scale)
#--> A3: Suggest width 212, height 336 pixels
#print "A4: Suggest width %i, height %i pixels" % (210 * mm / h_scale, 297 * mm / v_scale) 
#--> A4: Suggest width 150, height 237 pixels
def main(argv):
    if len(argv) != 3:
        sys.stderr.write('usage: dither_rbg.py organism_image single_chromosome.fasta\n')
        sys.exit()
    image = argv[1]
    seq_file = argv[2]
    if ' ' in image:
        sys.stderr.write('image name cannot contain spaces, please rename the file\n')
        sys.exit()
    if not os.path.isfile(image):
        sys.stderr.write('cannot find %s\n'%(image))
        sys.exit()
    if not os.path.isfile(seq_file):
        sys.stderr.write('cannot find %s\n'%(seq_file))
        sys.exit()
    
    pdf_file = os.path.splitext(os.path.basename(image))[0]+"_%s.pdf"

    if seq_file.endswith('.gz'):
        seq_fh = gzip.open(seq_file)
    else:
        seq_fh = open(seq_file)

    records = list(SeqIO.parse(seq_fh, "fasta"))
    records.sort(cmp=lambda x,y: cmp(len(y),len(x)))
    seq = str(records[0].seq)
    print len(seq)
    
    #Reduce runs of N to a single N to avoid visual distraction
    reduced = []
    reduced.append(seq[0])
    for i in seq[1:]:
        if i == 'N' and i == reduced[-1]:
            continue
        else:
            reduced.append(i)
    seq = Seq(''.join(reduced))

    print "Drawing %s using %s" % (image, seq_file)
    for name, shape, png_file in [
            ("A4", (150, 237), image),
            ("A3", (212, 336), image),
            ("A2", (300, 475), image),
            ("A1", (424, 672), image),
            ("A0", (600, 951), image),
        ]:
        if os.path.isfile(pdf_file % name):
            print "Skipping as %s exists..." % (pdf_file % name)
            continue
        print "Size %s, using %i by %i pixels from %s" \
            % (name, shape[0], shape[1], png_file)
        im = Image.open(png_file).resize(shape)
        run(im, seq, pdf_file % name, name)

if __name__ == '__main__':
    main(sys.argv)
    sys.exit()