#!/usr/bin/env python3

import spectral
import os
import matplotlib.pyplot as plt

bdr = spectral.open_image("data/outputs/prosail-bdr.bsq.hdr")

rgb = (260, 140, 20)

# This doesn't seem to work correctly...?
# Also, image might be flipped? May need to transform the R output.
spectral.imshow(bdr, rgb)
