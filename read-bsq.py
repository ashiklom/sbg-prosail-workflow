#!/usr/bin/env python3

import spectral
import os
import matplotlib.pyplot as plt

bdr = spectral.open_image("data/outputs/prosail-bdr.bsq.hdr")

rgb = (260, 140, 20)

spectral.imshow(bdr, rgb)
