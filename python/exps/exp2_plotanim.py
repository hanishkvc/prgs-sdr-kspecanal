#!/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

fig = plt.figure()

def update(frame):
    print("In Update")
    return plt.plot(np.random.rand(10))

def test_anim_blocking():
    print("Starting Anim")
    ani = FuncAnimation(fig, update, blit=True)
    plt.plot(np.random.rand(10))
    plt.show()

def test_anim_nonblocking(bWorkingLogic):
    print("Starting Anim")
    ani = FuncAnimation(fig, update, blit=True)
    plt.plot(np.random.rand(10))
    plt.show(block=False)
    #time.sleep(10)
    while True:
        print("In Main While")
        plt.plot(np.random.rand(10))
        #plt.draw()
        if (bWorkingLogic):
            plt.pause(1)
        else:
            time.sleep(5)

def test_nonblocking():
    plt.plot(np.random.rand(10))
    plt.show(block=False)
    while True:
        print("In Main While")
        plt.plot(np.random.rand(10))
        plt.pause(1)
        #time.sleep(5)

#test_anim_blocking()
test_anim_nonblocking(False)
#test_anim_nonblocking(True)
test_nonblocking()

