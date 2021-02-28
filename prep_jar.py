#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

global inv_kBT
inv_kBT = 1./0.5961610799306883

n_max=15

np.seterr(all='raise')

def smooth_f( frc, window  ):
    smooth_f      = np.zeros((len(frc),1))
    #smooth_f[:,0] = np.copy(frc[:])
    for i in range(len(frc)):
          w = window
          if i < w:
             w = i
          if i + w + 1>= len(frc):
             w = len(frc)-i-2
          smooth_f[i,0] = np.sum(frc[i-w:i+w+1])/(2*w+1)
    return smooth_f


def dWdx( x, ene ):
    smooth_f      = np.zeros((len(ene),1))
    if (x[1]-x[0]) != 0.:
        smooth_f[0,0] = (ene[1]-ene[0])/(x[1]-x[0])
    for i in range(1,len(ene)-1):
        if (x[i+1]-x[i-1]) != 0.:
          smooth_f[i,0] = (ene[i+1]-ene[i-1])/(x[i+1]-x[i-1])
    if (x[-1]-x[len(x)-2]) != 0.:
        smooth_f[-1,0] = (ene[-1]-ene[len(ene)-2])/(x[-1]-x[len(x)-2])
    return smooth_f

def exp_average( list_run_names ):
    count = 0
    minU  = 0.
    maxU  = 0.
    total = 0.
    for f in list_run_names:
        
        a = np.loadtxt(f)
        
        m = np.min(a[:,3])
        M = np.max(a[:,3])

        if m < minU:
            minU = m
        if M > maxU:
            maxU = M
        count += 1
        total += np.sum(a[:,3])/len(a[:,3])

    mean = total*inv_kBT/count
    midU = 0.5*(minU+maxU)*inv_kBT
    print("range of U is %e to %e, midpoint in kT: %e mean: %e\n" %\
          (minU, maxU, midU, mean))

    offSet = midU


    count = 0
    for f in list_run_names:
        print(f)

        a = np.loadtxt(f)

        a[:,3] *= inv_kBT
        a[:,3] -= offSet

        if count == 0:
            x = np.zeros((len(a),1))
            total        = np.zeros((len(a),1))
            try:
                total[:,0]  += np.exp(a[:,3])
            except Exception as e:
                print(e)
                exit(1)
                    
            ##save force versus extension:
            ##but actually, invisible string is getting shorter.
            x[:,0] = a[::-1,0]
        else:
            try:
                total[i,0]  += np.exp(a[:,3])
            except Exception as e:
                print(e)
                exit(1)
        count += 1.

    if count == 0:
        print( "error count == 0"+str(have))

    mean = (total / count) 

    fine = True
    try:
        mean  = (np.log(mean)+offSet) / inv_kBT
    except Exception as e:
        print(e)
        fine = False

    if not fine:    
       for i in range(len(mean)):
          print("%i %e" % (i,mean[i]))
          try:
             m  = np.log(mean[i])
             m += offSet
             m /= inv_kBT
          except Exception as e:
             print(e)
             print("%i %e" % (i,mean[i]))
             exit(1)
    return x, mean




if len(sys.argv) > 1:
    stem = sys.argv[1]
else:
    print("require a stem name")
    quit()

have=[]

for i in range(1,99):
    fName = stem+("/dist_vs_t_stretch%i_0.dat" % i)
    
    try:
       f     = open(fName, "r")
    except:
       continue
    print("found %s" % fName)

    joined_fName = stem + "/dist_vs_t_stretch%i_joined.dat" % i
    joined_f = open(joined_fName, "w")

    base_U = 0.
    latest = 0
    for j in range(1,n_max+1):
        fName = stem + "/dist_vs_t_stretch%i_%i.dat" % (i,j)
        try:
           f_array = np.loadtxt(fName)
           print("loaded: "+fName)
        except:
           print("Did not find: %s" % fName)
           continue
        
        f_array[:,3] += base_U
        base_U        = f_array[-1,3]

        if j == 1:
            joined = np.copy(f_array)
        else:
            joined = np.vstack( (joined, f_array) )

        latest = j          
               
        if latest >= n_max:
            have.append(joined_fName)
            ##save in append mode?
            smoothed_f = smooth_f( joined[:,2], 1000 )
            print(np.shape(joined))
            print(np.shape(smoothed_f))
            joined     = np.hstack( (joined, smoothed_f) )
            np.savetxt( joined_fName, joined )

            print( ("have %i runs for: " % latest) +str(have))
            break #enforce all runs same length

##make a file of distance_long X runs_wide work values 
if len(have) > 0:
    
    print("joining: "+str(have))

    count = 0
    for f in have:
        a = np.loadtxt(f)
        if count == 0:
            joined = np.zeros((len(a[:,3]),len(have)))
            joined[:,0] = np.copy(a[:,3])
        else:
            print(np.shape(joined))
            print(np.shape(a))
            joined[:,count] = np.copy(a[:,3])
            print("now:"+str(np.shape(joined)))
        count += 1
    np.savetxt("%s_allWorks.dat" % stem, joined)

