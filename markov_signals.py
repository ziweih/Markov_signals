from __future__ import division
import mdtraj as md
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import numpy as np
# from numpy import zeros, sqrt, square, log, where, pi, mean, average, arange, sin, cos
from math import sqrt, log, pi, sin, cos
from textwrap import wrap
from cycler import cycler
from glob import glob
import sys
from itertools import groupby
from decimal import Decimal
from mpmath import mp
import random


# Randomly generate integer of either 0 or 1
def random_bool():
    # R = random.randint(0,5,)
    R = random.randint(0,1)
    return R


# Randomly generate number between [0,1]
def random_range():
    R = random.uniform(0,1)
    return R


# Box-Muller transformation - generates independent gaussian numbers z1, z2
def gaussian():
    # z1 = sqrt(- 2 * log(random_seed())) * cos(2 * pi * random_seed())
    z1 = sqrt(- 2 * log(random_range())) * cos(2 * pi * random_range())
    #z2 = sqrt(- 2 * log(random()) * sin(2 * pi * random())
    #print iseed
    #print random(), random()
    # print 'z1 = ', z1
    return z1


# Generate a main Markovian signal
def main_signal():
    dt, n_decimals, t_f, n_frames, A, nstep, k_1, k_2, K = parameters()
    x = A * random_bool()        # starting position
    t = 0.0             # starting time
    
    main_signal_time = np.array([t])
    main_signal_traj = np.array([x])
    print 'MAIN SIGNAL \n'

    for i in xrange(nstep+1):
        if t <= t_f:
            R = random_range()
            print 'STEP ', i
            print 'x = ', x, '\nt = ', t
            print 'R = ', R

            if x == 0:
                t_B = - (1 / k_2) * log(R)
                print 't_B unrounded = ', t_B
                t_B = float(str(round(t_B, n_decimals)))
                num_data = int(round(t_B / dt))
                t_interval = np.linspace(dt, t_B, num=num_data)
                main_signal_traj = np.append(main_signal_traj, np.full(len(t_interval), x))
                main_signal_time = np.append(main_signal_time, t + t_interval)
                print 'Survival time of B = ', t_B
                print 'num_data = ', num_data
                print 't_interval = ', t_interval
                print 'appended times = ', t + t_interval
                print 'appended traj = ', np.full(len(t_interval), x)
                # print 'time array = ', main_signal_time
                # print 'traj array = ', main_signal_traj
                print 'length of t_interval =', len(t_interval)
                print 'length of traj,time = ', len(main_signal_time)
                x = x + A * random_bool()
                t = main_signal_time[-1]
            else:
                t_A = - (1 / k_1) * log(R)
                print 't_A unrounded = ', t_A
                t_A = float(str(round(t_A, n_decimals)))
                num_data = t_A / dt
                t_interval = np.linspace(dt, t_A, num=num_data)
                main_signal_traj = np.append(main_signal_traj, np.full(len(t_interval), x))
                main_signal_time = np.append(main_signal_time, t + t_interval)
                print 'Survival time of A = ', t_A
                print 'num_data = ', num_data
                print 't_interval = ', t_interval
                print 'appended times = ', t + t_interval
                print 'appended traj = ', np.full(len(t_interval), x)
                # print 'time array = ', main_signal_time
                # print 'traj array = ', main_signal_traj
                print 'length of t_interval =', len(t_interval)
                print 'length of traj and time = ', len(main_signal_time)
                x = x - A * random_bool()
                t = main_signal_time[-1]
        else:
            break
        print
    main_signal_time = main_signal_time[0:n_frames+1]
    main_signal_traj = main_signal_traj[0:n_frames+1]
    np.save('data/main_signal.npy', zip(main_signal_time, main_signal_traj))
    print 'Summary:'
    print 'time = ', main_signal_time
    print 'traj = ', main_signal_traj
    print 'last time step in saved file = ', main_signal_time[-1]
    print 'length of saved traj and time', len(main_signal_time)
    print
    print 't_f = ', t_f
    print 'n_frames = ', n_frames
    print 'n_decimals = ', n_decimals
    print
    plt.figure(figsize=(20,10))
    plt.plot(main_signal_time[0:n_frames+1], main_signal_traj[0:n_frames+1])
    plt.xlim(0,t_f)
    plt.savefig('figures/main_signal.pdf')
    # plt.show()
    return main_signal_time, main_signal_traj


# Generate unimportant signals with shorter survival times than the main signal
    # and with amplitude lower or equal to that of main signal
def other_signal():
    num_signals = 5
    other_signal_traj = {}
    other_signal_time = {}
    
    for i in xrange(num_signals):
        dt, n_decimals, t_f, n_frames, A, nstep, k_1, k_2, K = parameters()
        A = random.uniform(0, A/2)
        k_1 = random.uniform(k_1*10,100)
        k_2 = random.uniform(k_2*10,100)
        x = A * random_bool()               # starting position
        t = 0.0
        print 'OTHER SIGNAL ', i
        print 'parameters: A = ', A, ', k_1 = ', k_1, ', k_2 = ', k_2
        print 'initial: x = ', x, ', t = ', t, ', dt = ', dt

        time_temp = np.array([t])
        traj_temp = np.array([x])
        print 'initial time and traj list = ', time_temp, traj_temp

        for j in xrange(nstep+1):
            if t <= t_f:
                R = random_range()
                print 'STEP ', j
                print 't = ', t, '\nx = ', x
                print 'R = ', R

                if t != time_temp[-1]:
                    print 'ERROR FOR SIGNAL ', i, 'AT STEP ', j

                if x == 0:
                    t_B = - (1 / k_2) * log(R)
                    print 't_B unrounded = ', t_B
                    t_B = float(str(round(t_B, n_decimals)))
                    num_data = t_B / dt
                    t_interval = np.linspace(dt, t_B, num=num_data)
                    traj_temp = np.append(traj_temp, np.full(len(t_interval), x))
                    time_temp = np.append(time_temp, t + t_interval)
                    print 'Survival time of B = ', t_B
                    print 'num_data = ', num_data
                    print 't_interval = ', t_interval
                    print 'appended times = ', t + t_interval
                    print 'appended traj = ', np.full(len(t_interval), x)
                    # print 'time array = ', time_temp
                    # print 'traj array = ', traj_temp 
                    print 'length of t_interval =', len(t_interval)
                    print 'length of traj and time = ', len(time_temp)
                    x = x + A * random_bool()
                    t = time_temp[-1]
                else:
                    t_A = - (1 / k_1) * log(R)
                    print 't_A unrounded = ', t_A
                    t_A = float(str(round(t_A, n_decimals)))
                    num_data = t_A / dt
                    t_interval = np.linspace(dt, t_A, num=num_data)
                    traj_temp = np.append(traj_temp, np.full(len(t_interval), x))
                    time_temp = np.append(time_temp, t + t_interval)
                    print 'Survival time of A = ', t_A
                    print 'num_data = ', num_data
                    print 't_interval = ', t_interval
                    print 'appended times = ', t + t_interval
                    print 'appended traj = ', np.full(len(t_interval), x)
                    # print 'time array = ', time_temp
                    # print 'traj array = ', traj_temp 
                    print 'length of t_interval =', len(t_interval)
                    print 'length of traj and time = ', len(time_temp)
                    x = x - A * random_bool()
                    t = time_temp[-1]
            else:
                break
            print
        print
        print 'Summary'
        print 'uncut length of each saved traj and time = ', len(time_temp), len(traj_temp)
        time_temp = time_temp[0:n_frames+1]
        traj_temp = traj_temp[0:n_frames+1]
        np.save('data/other_signal_{}.npy'.format(i), zip(time_temp,traj_temp))
        other_signal_time['{}'.format(i)] = time_temp
        other_signal_traj['{}'.format(i)] = traj_temp
        print 'time = ', time_temp
        print 'traj = ', traj_temp
        print 'last time step in saved file = ', time_temp[-1]
        print 'length of each saved traj and time', len(time_temp), len(other_signal_time['{}'.format(i)])
        print
    print 't_f = ', t_f
    print 'n_frames = ', n_frames
    print 'n_decimals = ', n_decimals
    print
    plt.figure(figsize=(20,10))
    for i in xrange(num_signals):
        print 'lengths check: ', len(other_signal_time['{}'.format(i)])
        plt.plot(other_signal_time['{}'.format(i)][0:n_frames+1],\
                 other_signal_traj['{}'.format(i)][0:n_frames+1])
    plt.savefig('figures/other_signal.pdf')
    # plt.show()
    return other_signal_time, other_signal_traj


# Generate random Gaussian noise
def noise():
    num_signals = 10
    noise_traj = {}
    noise_time = {}
    
    for i in xrange(num_signals):
        g = gaussian()
        dt, n_decimals, t_f, n_frames, A, nstep, k_1, k_2, K = parameters()
        x = g      # starting position
        t = 0.0
        print 'RANDOM NOISE ', i
        print 'parameters: A = ', A, 'initial g = ', g, ', k_1 = ', k_1, ', k_2 = ', k_2
        print 'initial: x = ', x, ', t = ', t, ', dt = ', dt

        time_temp = np.array([t])
        traj_temp = np.array([x])

        for j in xrange(nstep+1):
            if t <= t_f:
                g = gaussian()
                print 'STEP ', j
                print 't = ', t, '\nx = ', x
                print 'g = ', g
                t = t + dt
                x = g
                time_temp = np.append(time_temp, t)
                traj_temp = np.append(traj_temp, x)
                print 'new x = ', x, '\nnew t = ', t
                print 'length of traj and time = ', len(time_temp)
            else:
                break
            print
        print
        time_temp = time_temp[0:n_frames+1]
        traj_temp = traj_temp[0:n_frames+1]
        np.save('data/noise_{}.npy'.format(i), zip(time_temp,traj_temp))
        noise_time['{}'.format(i)] = time_temp
        noise_traj['{}'.format(i)] = traj_temp
        print 'Summary:'
        print 'time = ', time_temp
        print 'traj = ', traj_temp
        print 'last time step in saved file = ', time_temp[-1]
        print 'length of each saved traj and time', len(time_temp)
        print
    print 't_f = ', t_f
    print 'n_frames = ', n_frames
    print 'n_decimals = ', n_decimals
    print

    plt.figure(figsize=(20,10))
    for i in xrange(num_signals):
        plt.plot(noise_time['{}'.format(i)][0:n_frames+1],\
                 noise_traj['{}'.format(i)][0:n_frames+1])
    plt.savefig('figures/noise.pdf')
    plt.show()
    return noise_time, noise_traj
    

# Parameters for main signal
def parameters():
    dt = 0.001
    n_decimals = 3             # number of decimal places in dt
    t_f = 10.0                # final time of signal (want all traj to be same length)
    n_frames = int(t_f / dt)   # number of frames to write in signal trajectories
    A = 5.0                    # amplitude of signal
    nstep = 10000000           # number of steps
    k_1 = 4                    # forward rate of conversion A -> B
    k_2 = 4                    # backward rate              A <- B
    K = k_2 / k_1
    # print dt, n_decimals, t_f, n_frames, A, nstep, k_1, k_2, K
    return dt, n_decimals, t_f, n_frames, A, nstep, k_1, k_2, K


# random_bool()
# random_range()
# gaussian()
# main_signal()
other_signal()
# noise()
# parameters()






