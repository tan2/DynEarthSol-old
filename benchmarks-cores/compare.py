#!/usr/bin/python3
from __future__ import print_function
import sys, os
import numpy as np
sys.path.append('../')
from Dynearthsol import Dynearthsol

def first_invariant(t):
    nstr = t.shape[1]
    ndims = 2 if (nstr == 3) else 3
    return np.sum(t[:,:ndims], axis=1) / ndims


def second_invariant(t):
    '''The second invariant of the deviatoric part of a symmetric tensor t,
    where t[:,0:ndims] are the diagonal components;
      and t[:,ndims:] are the off-diagonal components.'''
    nstr = t.shape[1]

    # second invariant: sqrt(0.5 * t_ij**2)
    if nstr == 3:  # 2D
        return np.sqrt(0.25 * (t[:,0] - t[:,1])**2 + t[:,2]**2)
    else:  # 3D
        a = (t[:,0] + t[:,1] + t[:,2]) / 3
        return np.sqrt( 0.5 * ((t[:,0] - a)**2 + (t[:,1] - a)**2 + (t[:,2] - a)**2) +
                        t[:,3]**2 + t[:,4]**2 + t[:,5]**2)


class Stuff():
    pass


def read_data(des, frame):
    stuff = Stuff()

    stuff.T = des.read_field(frame,'temperature')
    coordinate = des.read_field(frame, 'coordinate')
    stuff.x = np.array(coordinate[:,0])
    stuff.z = np.array(coordinate[:,-1])
    velocity = des.read_field(frame, 'velocity')
    stuff.vx = np.array(velocity[:,0])
    stuff.vz = np.array(velocity[:,-1])
    stuff.pls = des.read_field(frame,'plastic strain')

    stress = des.read_field(frame, 'stress')
    stuff.tI = first_invariant(stress)
    stuff.tII = second_invariant(stress)
    strain = des.read_field(frame, 'strain')
    stuff.sI = first_invariant(strain)
    stuff.sII = second_invariant(strain)
    strain_rate = des.read_field(frame, 'strain-rate')
    stuff.srI = first_invariant(strain_rate)
    stuff.srII = second_invariant(strain_rate)

    stuff.visc = des.read_field(frame, 'viscosity')

    marker_data = des.read_markers(frame, markersetname)
    field = marker_data[markersetname + '.coord']
    stuff.m_x = field[:,0]
    stuff.m_z = field[:,1]
    stuff.m_id = marker_data[markersetname + '.id']
    stuff.m_mat = marker_data[markersetname + '.mattype']
    stuff.m_time = marker_data[markersetname + '.time']
    return stuff


def reldiff(oldf, newf):
    m = np.abs(oldf).max()
    diff = np.abs(newf - oldf)
    if m == 0.:
        return diff.max(), diff.std()
    else:
        return diff.max()/m, diff.std()/m


def show_msg(kind,max,sigma):
    if max+sigma > 1.e-8:
        print('  %s:\t\t%.3e %.3e (> 1.e-8)'%(kind,max, sigma))
        inc = 1
    else:
        print('  %s:\t\t%.3e %.3e'%(kind,max, sigma))
        inc = 0
    return inc


def compare(old, new):
    inc = 0

    max, sigma = reldiff(old.T, new.T)
    inc += show_msg('Temperature',max,sigma)

    max, sigma = reldiff(old.x, new.x)
    inc += show_msg('X coordinate',max,sigma)

    max, sigma = reldiff(old.z, new.z)
    inc += show_msg('Z coordinate',max,sigma)

    max, sigma = reldiff(old.vx, new.vx)
    inc += show_msg('X velocity',max,sigma)

    max, sigma = reldiff(old.vz, new.vz)
    inc += show_msg('Z velocity',max,sigma)

    max, sigma = reldiff(old.pls, new.pls)
    inc += show_msg('Pl. strain',max,sigma)

    max, sigma = reldiff(old.tI, new.tI)
    inc += show_msg('Stress I',max,sigma)

    max, sigma = reldiff(old.tII, new.tII)
    inc += show_msg('Stress II',max,sigma)

    max, sigma = reldiff(old.sI, new.sI)
    inc += show_msg('Strain I',max,sigma)

    max, sigma = reldiff(old.sII, new.sII)
    inc += show_msg('Strain II',max,sigma)

    max, sigma = reldiff(old.srI, new.srI)
    inc += show_msg('S. rate I',max,sigma)

    max, sigma = reldiff(old.srII, new.srII)
    inc += show_msg('S. rate II',max,sigma)

    max, sigma = reldiff(old.visc, new.visc)
    inc += show_msg('Viscosity',max,sigma)

    max, sigma = reldiff(old.m_x, new.m_x)
    inc += show_msg('Marker X',max,sigma)

    max, sigma = reldiff(old.m_z, new.m_z)
    inc += show_msg('Marker Z',max,sigma)

    max, sigma = reldiff(old.m_mat, new.m_mat)
    inc += show_msg('Marker Mat',max,sigma)

    max, sigma = reldiff(old.m_time, new.m_time)
    inc += show_msg('Marker Time',max,sigma)

    return inc


olddir = sys.argv[1]
curdir = os.getcwd()

if len(sys.argv) > 3:
    frame = int(sys.argv[3])
    newdir = sys.argv[2]
    modelname = 'result'
else:
    frame = int(sys.argv[2])
    newdir = curdir
    modelname = 'benchmark'

# name holder
old = 0
new = 0

markersetname = 'markerset'

try:
    # read old and new results

    os.chdir(olddir)
    des = Dynearthsol(modelname)
    old = read_data(des, frame)

    os.chdir(newdir)
    des = Dynearthsol(modelname)
    new = read_data(des, frame)

    # compare results
    print()
    print('Relative difference (max, stddev) of frame =', frame,
          ' step =', int(des.steps[frame]))
    print('  ---')
    inc = compare(old, new)
    print('')
    if inc == 0:
        print('  Status: Normal round-off error~')
    else:
        print('  Status: !!!!!!!!!! SOMETHING WRONG !!!!!!!!!!')
    print('  ---')


finally:
    # restort to original directory
    os.chdir(curdir)
