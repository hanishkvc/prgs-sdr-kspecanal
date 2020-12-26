#!/usr/bin/env python3

import matplotlib.pyplot as plt

def handle_clicked(event):
    print("Btn clicked...")

f1 = plt.figure(figsize=(6,4))
gs = f1.add_gridspec(2,2)
ax = f1.add_subplot(gs[0,0])
btn = plt.Button(ax, "test")
btn.on_clicked(handle_clicked)
plt.show()

