from __future__ import division
import h5py
import sys, imp
import numpy as np
from numpy import rec
import math
import matplotlib
from matplotlib import pyplot as plt
from scipy import stats
import scipy
import cPickle as pickle
import gizmo_analysis as gizmo
import utilities as ut
import sat_fct as sf
from all_fcts import *
import matplotlib
from collections import Counter
from astropy.cosmology import FlatLambdaCDM

h = 0.702
Fs = 10
fontsize = 15

# Nice plots - thick axes
plt.rc('axes',linewidth=1.5)

#Dictionary of chemical species ids
spec_dict = {'Z': 0, 'H': 0, 'He': 1, 'C': 2, 'N': 3, 'O': 4, 'Ne': 5, 
		'Mg': 6, 'Si': 7, 'S': 8, 'Ca': 9, 'Fe': 10}

#From Anders & Grevesse 1989, Fe from Sneden et al 1992
solar_abnds = {'Mg': 7.58, 'Ca': 6.36, 'Si': 7.55, 'Ti': 4.99, 'Fe': 7.52, 'O': 8.93,
			   'H': 12.}
			   
#Atomic mass units
amu = {'H': 1.008, 'He': 4.003, 'Mg': 24.305, 'Ca': 40.078, 'Si': 28.086, 'Ti': 47.88,\
	   'Fe': 55.933, 'O': 16.00}

# ## PARENT CLASS
###########################################################################################

class Particles:
    #initialize
    def __init__(self, workpath, fname, snapnum, part, header, data, runname, ldict, cdict, ind):
        self.fname = "%s/output/fname"%(workpath)
#         self.fname = "%s/fname"%(workpath) # GP
        self.snapnum = snapnum
        self.snap ='{:04d}'.format(snapnum)
        self.part = part
        self.header = header
        self.data = data
        self.runname = runname
        self.ind = ind     
        if ind == 0: # central/main
            self.color = cdict["%s_cen"%(runname)]
            self.label = ldict["%s_cen"%(runname)]
            self.symbol = '*'
            self.msize = 18
        elif ind == 1: 
        	self.color = cdict["%s_sat"%(runname)]
            self.label = ldict["%s_sat"%(runname)]
            self.symbol = 's'
            self.msize = 14
        elif ind == 2:
            self.color = 'salmon'
            self.label = None
            self.symbol = '<'
            self.msize = 14
        elif ind == 3:
            self.color = 'lightsalmon'
            self.label = None
            self.symbol = 'h'
            self.msize = 14
        elif ind == 4: ###### SATELLITE
        	self.color = 'g'
        	self.label = 'Satellite'
        	self.symbol = '^'
        	self.msize = 14
        
        self.mvir = self.data.field('Mvir(4)')[ind]/h
        self.rvir = self.data.field('Rvir(12)')[ind]/h
        self.mstar = self.data.field('M_star(65)')[ind]/h
        self.rmax = self.rvir*0.15
        
        self.centered = False
        self.n = 0 # will be set in child class
        self.pcen = part.host_positions[0]
        self.vcen = part.host_velocities[0]
        self.arrays = np.ndarray(self.n) # will be set in child class
        self.pos = [] # will be set in child class
        self.vel = [] # will be set in child class
        self.hh = self.header['hubble']
#         self.Hsol = 0.7381
        self.Zsol = 0.02
#         self.Fesol = 1.73e-3
#         self.alphasol = 2.22e-3 + 9.31e-4 + 1.08e-3 + 6.44e-4 + 1.01e-4 # Ne,Mg,Si,S,Ca
        self.booly = [] # initialize to all true in child classes, used to restrict particles

        self.get_arrays_from_file()
        
    def save(self, filename):
    	"""
    	
    	Uses pickle to save object
    	"""
        with open(filename, 'wb') as output:  # Overwrites any existing file.
            pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)

    def get_arrays_from_file(self):
        pass # these will be overwritten

    # Recenter can be passed r,m,v of all particles, or just a certain type
    # this will be decided in the execute file. Phase space arrays reset with new values in 
    # this function.
    def recenter(self, newcen = None, newvcen = None, fromscratch = True):
        """
        Recenter MUST take at least one new center, otherwise nothing done
        
        """
        
        if newcen is None:
        	newcen = self.part.host_positions[0]
        if newvcen is None:
        	newvcen = self.part.host_velocities[0]
        
        if self.centered == False:
			if fromscratch == True:
				self.reset_to_orig(recenter = True)

				self.pcen = np.array([0,0,0])
				self.pos = self.pos - newcen
				self.arrays.r = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
				self.arrays.x = self.pos[:,0]
				self.arrays.y = self.pos[:,1]
				self.arrays.z = self.pos[:,2]
			
				self.vcen = np.array([0,0,0])
				self.vel = self.vel - newvcen
				self.arrays.velx = self.vel[:,0]
				self.arrays.vely = self.vel[:,1]
				self.arrays.velz = self.vel[:,2]
			self.centered = True
    
    # This function puts all arrays in order of increasing radial distance (r)
    def orderp(self):
        pass # overwritten in child classes
    
    # Color set
    def set_clr(self, clr):
        self.color = clr
    
    # This method plots density profile. 
    def plot_rho(self, nbins=10, logbins = False, logR=False, clr = 'k', lw = 2, ls = '-', labels = None):
        PI = 3.14159265359
        self.orderp() # Get everything in order of increasing r
        rad = self.arrays.r[self.booly]
#         print np.amax(rad)
        mpart = self.arrays.mass[self.booly]
        rbins = gen_bins(rad, nbins, logbins)
    
        vperbin = np.ndarray(nbins) # will hold volume per bin
        mperbin = np.ndarray(nbins) # will hold mass per bin
        nperbin = np.ndarray(nbins) # will hold number of particles per bin
        midbin = np.ndarray(nbins) # will store middle of r bin
    
        # Calculate volume per bin
        for c in range(nbins):
            vperbin[c] = (rbins[c+1] ** 3 - rbins[c] ** 3) * (4./3) * PI # volume per bin (kpc^3)
            radinbin = rad[((rad >= rbins[c]) & (rad < rbins[c+1]))]
            massinbin = mpart[((rad >= rbins[c]) & (rad < rbins[c+1]))] # nec for diff masses
            nperbin[c] = np.size(radinbin) # number of particles in each bin
            mperbin[c] = np.sum(massinbin)
            midbin[c] = (rbins[c+1] + rbins[c]) / 2

        # Volume density = npart * mpart / vbin
        rho = mperbin / vperbin    
        if logR == True:
        	plt.plot(midbin, rho, color = clr, linewidth = lw, linestyle = ls, label = labels) # Pass string for legend
        	plt.xlabel(r'$\rm R (kpc)$', fontsize = 22)
        	plt.loglog()
        else:
        	plt.plot(midbin, np.log10(rho), color = clr, linewidth = lw, linestyle = ls, label = labels) # Pass string for legend
        	plt.xlabel(r'$\rm R (kpc)$', fontsize = 22)
        if labels is not None:
        	plt.legend(loc = 1, numpoints = 1, frameon = False, fontsize = 22)
        plt.ylabel(r'$\rm \rho(R) (M_{\odot}/kpc^{3})$', fontsize = 22)
        return midbin,rho
        
        
        
    # This method plots mass profile
    def m_profile(self, nbins=100, cumu = True, logbins = False, logR=False, clr = None):
        self.orderp() # Get everything in order of increasing r
        rad = self.arrays.r[self.booly]
#         print np.amax(rad)
        mpart = self.arrays.mass[self.booly]
        labels = self.label
        rbins = gen_bins(rad, nbins, logbins)

        mperbin = np.ndarray(nbins) # will hold mass per bin
        nperbin = np.ndarray(nbins) # will hold number of particles per bin
        midbin = np.ndarray(nbins) # will store middle of r bin
        
        if cumu == False:
            # Calculate mass per bin
            for c in range(nbins):
                radinbin = rad[((rad >= rbins[c]) & (rad < rbins[c+1]))]
                massinbin = mpart[((rad >= rbins[c]) & (rad < rbins[c+1]))] # nec for diff masses
                nperbin[c] = np.size(radinbin) # number of particles in each bin
                mperbin[c] = np.sum(massinbin)
                midbin[c] = (rbins[c+1] + rbins[c]) / 2
        else:
            # Calculate cumu mass per bin
            massinbin = mpart[((rad >= rbins[0]) & (rad < rbins[1]))]
            mperbin[0] = np.sum(massinbin)
            midbin[0] = (rbins[1] + rbins[0]) / 2
            for c in range(nbins-1):
                massinbin = mpart[((rad >= rbins[c+1]) & (rad < rbins[c+2]))] # nec for diff masses
                mperbin[c+1] = np.sum(massinbin) + mperbin[c]
                midbin[c+1] = (rbins[c+2] + rbins[c+1]) / 2
#             print np.log10(amax(mperbin))
            
        plt.figure(1, figsize=(7, 6))
        if logR == True:
            plt.plot(np.log10(midbin), np.log10(mperbin), color = clr, linewidth = 2, label = labels) # Pass string for legend
            plt.xlabel(r'$\rm log \,R (kpc)$', fontsize = 22)
        else:
            if clr is None:
                plt.plot(midbin, np.log10(mperbin), color = self.color, linewidth = 2, label = labels) # Pass string for legend
            else:
            	plt.plot(midbin, np.log10(mperbin), color = clr, linewidth = 2, label = labels) # Pass string for legend
            plt.xlabel(r'$\rm R (kpc)$', fontsize = 22)
            plt.legend(loc = 4, numpoints = 1, frameon = False, fontsize = 12)
            plt.ylabel(r'$\rm log \,M (M_{\odot})$', fontsize = 22)
            plt.xticks(weight = 'bold')
            plt.yticks(weight = 'bold')
		
    def half_mass_r(self, r = 5.0, radd = 0.5):
        """
        Radius within which radd fraction of mass is enclosed
        returns rhalf, defaults to radd = 50% (half mass radius)
	
        """
	   	# First get all particle arrays within stated distance
        self.recenter(newcen = self.part.host_positions[0], newvcen = self.part.host_velocities[0])
        mass = self.arrays.mass[self.arrays.r < r]
        r = self.arrays.r[self.arrays.r < r]
        masstot = np.sum(mass)
	
        masssort = mass[np.argsort(r)]
        rr = r[np.argsort(r)]
        N = len(mass)
        massarray = np.zeros(N)
        massarray[0] = masssort[0]
        for t in range(N-1):
            massarray[t+1] = (massarray[t] + masssort[t+1])
        massarr = massarray  / masstot # Divide out total mass to get mass fraction
        # Calculate half mass radius
        rhalf = rr[np.argmin(np.abs(massarr - radd))]
        return rhalf

    def half_mass_r_proj(self, r, radd = 0.5, face = 'xy'):
        """
        2D radius within which radd fraction of mass is enclosed
        returns rhalf, defaults to radd = 50% (half mass radius)
	
        """
	   	# First get all particle arrays within stated distance
        if face == 'xy':
            r2d = np.sqrt(self.arrays.x**2 + self.arrays.y**2)
        elif face == 'xz':
            r2d = np.sqrt(self.arrays.x**2 + self.arrays.z**2)
        else:
            r2d = np.sqrt(self.arrays.y**2 + self.arrays.z**2)
        mass = self.arrays.mass[r2d < r]
        r2d = r2d[r2d < r]
        masstot = np.sum(mass)
	
        masssort = mass[np.argsort(r2d)]
        rr = r2d[np.argsort(r2d)]
        N = len(mass)
        massarray = np.zeros(N)
        massarray[0] = masssort[0]
        for t in range(N-1):
            massarray[t+1] = (massarray[t] + masssort[t+1])
        massarr = massarray  / masstot # Divide out total mass to get mass fraction
        # Calculate half mass radius
        rhalf2d = rr[np.argmin(np.abs(massarr - radd))]
        return rhalf2d

    def mass_in_r(self, r = None):
        """
        mass within r
	
        """
        self.recenter(newcen = self.part.host_positions[0], newvcen = self.part.host_velocities[0])
        
        # first restrict to rvir if no other r given
        
        if r is not None:
            rad = self.arrays.r[self.arrays.r < r]
            mass = self.arrays.mass[self.arrays.r < r]
        else:
            rad = self.arrays.r[self.arrays.r < self.rvir]
            mass = self.arrays.mass[self.arrays.r < self.rvir]
		
        # Now order for increasing r
        msort = mass[np.argsort(rad)]
        rsort = rad[np.argsort(rad)]
        
        m_in_r = np.array([np.sum(msort[0:i]) for i in range(len(rsort))])
			
        return rsort, m_in_r
        
    def v_with_r(self, r = None, h = 2):
        """
        velocity with r
        CHECK
	
        """
        self.recenter(newcen = self.part.host_positions[0], newvcen = self.part.host_velocities[0])
        
        # first restrict to close to disk

        rad = self.arrays.x[np.abs(self.arrays.y) < h]
        vel = self.arrays.velz[np.abs(self.arrays.y) < h]
        
        if r is not None:
            vel = vel[((rad > 0)&(rad < r))]
            rad = rad[((rad > 0)&(rad < r))]
        else:
            vel = vel[((rad > 0)&(rad < self.rvir))]
            rad = rad[((rad > 0)&(rad < self.rvir))]

		
        # Now order for increasing r
        vsort = vel[np.argsort(rad)]
        rsort = rad[np.argsort(rad)]
        
			
        return rsort, vsort

    def plot_mass_2Dhist(self, r_ext = 50, h_ext = 50, bins = 500, cmap = 'bone', cmin = 0, topdown = True):
        """
        This function plots a visualization of the particle type from the z axis looking down,
            weighted by the particle mass
        
        
        """
        fig = plt.figure(1, figsize = (10,10))
        ax1 = fig.add_subplot(1,1,1)
        ax1.patch.set_facecolor('black')
        
        if topdown == True:
        	rr = np.sqrt((self.arrays.x[self.booly] **2) + (self.arrays.y[self.booly] **2))
        else:
        	rr = np.sqrt((self.arrays.x[self.booly] **2) + (self.arrays.z[self.booly] **2))

        
        # Only choose those particles in selected region
        x = self.arrays.x[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]
        y = self.arrays.y[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]
        z = self.arrays.z[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]
        weights = self.arrays.mass[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]
        
        plot_m_2Dhist(x,y,z,weights,bins,cmap,cmin,topdown)
        
    def plot_mass_2Dhist_basic(self, bins = 500, cmap = 'bone', cmin = 0, topdown = True, xylims = None):
        """
        This function plots a visualization of the particle type from the z axis looking down,
        weighted by the particle mass
        
        
        """
        x = self.arrays.x[self.booly]
        y = self.arrays.y[self.booly]
        z = self.arrays.z[self.booly]
        weights = self.arrays.mass[self.booly]
        if topdown == True:
            plt.hist2d(x, y, norm=matplotlib.colors.LogNorm(), bins = bins, weights = weights, cmap = cmap, cmin = cmin)
        else:
            plt.hist2d(x, z, norm=matplotlib.colors.LogNorm(), bins = bins, weights = weights, cmap = cmap, cmin = cmin)
        if xylims is not None:
            plt.xlim(xylims[0],xylims[1])
            plt.ylim(xylims[2],xylims[3])
        plt.show()
	
	
	
        
    def get_order(self):
        """
        Returns ordered indices of r
	
        """    
        return np.argsort(self.arrays.r) # array indices in order of increasing r
        
    def restrict_particles(self, r_restr_lo, r_restr_hi):
        """
        This function limits all particles to the specified annulus
        For the height, option is given to get separate regions above and below the disk(habs == True)
        or to simply go between two heights
    
        """
    
        # Only choose those particles in selected region

        self.booly = ((self.arrays.r > r_restr_lo) & (self.arrays.r < r_restr_hi))

    def restrict_parts(self, r_restr):
        """
        This function limits all particles to the specified sphere
    
        """
    
        # Only choose those particles in selected region
        self.booly = (self.arrays.r < r_restr)
        
    def reset_to_orig(self, recenter = False):
        """
        This function reloads the whole file using the snapnum
        
        """
        self.booly=np.zeros(len(self.arrays.id)) == 0
        if recenter == True:
            self.get_arrays_from_file() 
            self.pcen = self.part.host_positions[0]
            self.vcen = self.part.host_velocities[0]
        self.centered = False
        
        
        
    def plot_aplha_Fe_vs_Fe_H_hist(self, binx = 300, biny = 300, lims = [-1,-1]):
        pass
    
    def return_non_dup_ids(self):
        """
        This function removes all duplicate ids and returns the list of ids without the duplicates
        
        """
        new_ids = np.array([k for k,v in Counter(self.arrays.id).iteritems() if v == 1])
        return new_ids
    
    def get_r_from_ids(self, idlist):
        """
        This function returns the current radius of all particles with ids in idlist
        If a certain id is not found, returns -1 (then a later snapshot must be searched outside
            this function)
    
        """
        sorted_ids = self.arrays.id[argsort(self.arrays.id)]
        r_now = np.ndarray(len(idlist))
        for t in range(len(idlist)):
            r_now[t] = self.arrays.r[np.where(self.arrays.id == new_ids[t])][0]
        return r_now

    def calc_mass_within_r(self, rad, recenter = False):
        """
        This function returns total mass of given species within specified radius
        """
        if self.centered == False:
        	self.recenter(newcen = self.part.host_positions[0], newvcen = self.part.host_velocities[0])
        return np.sum(self.arrays.mass[self.booly][self.arrays.r[self.booly] < rad])        
        
    def plot_cumu_mass(self, rad = 5, normed = False):
        """

        """
	   	# First get all particle arrays within stated distance
        self.recenter(newcen = self.part.host_positions[0], newvcen = self.part.host_velocities[0])
        mass = self.arrays.mass[self.arrays.r < rad]
        r = self.arrays.r[self.arrays.r < rad]
        masstot = np.sum(mass)
	
        masssort = mass[np.argsort(r)]
        N = len(mass)
        massarray = np.zeros(N)
        massarray[0] = masssort[0]
        for t in range(N-1):
            massarray[t+1] = (massarray[t] + masssort[t+1])
        massarr = massarray  / masstot # Divide out total mass to get mass fraction
        
        if normed == True:
        	return r,massarr
        else:
        	return r,massarray

    def calc_sig(self, r = None, ttype = 'los', mthd = 'iqr', axis = 'x'):
        """
        Calculates velocity dispersion
        r = restrict to this radius
        ttype = line of sight or 3d (if 3d ignores axis)
        mthd - either straight sigma or interquartile range (more robust to outliers)
        axis - for los method, which axis is along los
	
	
        """
        if r is not None: # must restrict particles to this radius
            self.restrict_parts(r)
        if ttype == 'los':
            if axis == 'x':
                vel = self.arrays.velx[self.booly]
            elif axis == 'y':
                vel = self.arrays.vely[self.booly]
            else:
                vel = self.arrays.velz[self.booly]
        else: # 3d velocity
            vel = np.sqrt(self.arrays.velx[self.booly] **2 + self.arrays.vely[self.booly] **2 + self.arrays.velz[self.booly] **2)
            
        if mthd == 'sig':
            return np.std(vel)
        else:
            return np.iqr(vel)
    
# ## CHILD CLASSES
###########################################################################################
class DarkParticles(Particles):
    # As a part of the initialization, get all arrays
    def get_arrays_from_file(self):
        idhalo=self.part['dark']['id'][:] # the ids of the halo particles
        mhalo=self.part['dark']['mass'][:]
        self.pos=self.part['dark']['position'][:]
        self.vel=self.part['dark']['velocity'][:]
        rhalo = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.n=len(mhalo) # number of halo particles
        self.arrays = rec.array([idhalo, mhalo, self.pos[:,0], self.pos[:,1], self.pos[:,2], rhalo, self.vel[:,0], \
            self.vel[:,1], self.vel[:,2]], dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'), ('r', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8')])
        self.booly=np.zeros(len(self.arrays.id)) == 0 # initialize to all true, used to restrict particles
    # Put all arrays in order of increasing distance from center (argsort r)    
    def orderp(self):
        order = np.argsort(self.arrays.r)
        self.arrays.r = self.arrays.r[order]
        self.arrays.id = self.arrays.id[order]
        self.arrays.mass = self.arrays.mass[order]
        self.arrays.x = self.arrays.x[order]
        self.arrays.y = self.arrays.y[order]
        self.arrays.z = self.arrays.z[order]
        self.arrays.velx = self.arrays.velx[order]
        self.arrays.vely = self.arrays.vely[order]
        self.arrays.velz = self.arrays.velz[order]
        self.booly = self.booly[order]
        
    def plot_Fe_H_vs_aplha_Fe_hist(self):
        print "ERROR: No Metallicity for Dark Matter"
        
    # Label for Dark Matter
    def get_labels(self):
        return 'Dark Matter'

    def plot_rho_power(self, r_restr = 500, nbins = 100, alpha = 0.2):
        """
        
        """
        self.recenter()
        self.restrict_parts(r_restr)
        r = self.arrays.r[self.booly]
        rsort = r[np.argsort(r)]
        power = rsort[199]
        rad,rho = self.plot_rho(nbins = nbins, logbins = True, logR = True, clr = self.color, lw = 2.5, labels = self.label)
        plt.axvspan(.001, power, facecolor=self.color, edgecolor = self.color, alpha=alpha)
        return rad,rho

    def downsize(self):
        """
        Completely cuts size down to restricted (by booly) to save space
        
        """
        self.pos=self.pos[self.booly]
        self.vel=self.vel[self.booly]
        self.arrays = self.arrays[self.booly]
        self.booly = self.booly[self.booly]
        self.n=len(self.booly) # number of halo particles
		

#############################################        
class StarParticles(Particles):
    def get_arrays_from_file(self):
        age = self.part['star']['form.time'][:] # cosmological
        lookback_time = 13.8 - age
        a = self.part['star'].prop('form.scalefactor')
#         print "Age is formation time"    
        idstar=self.part['star']['id'][:] # the ids of the star particles
        mstar=self.part['star']['mass'][:]
        self.pos=self.part['star']['position'][:]
        self.vel=self.part['star']['velocity'][:]
        rstar = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.Zstar=self.part['star']['massfraction'][:]
        self.Fe_H = get_abnd_ratios(self.Zstar, ('Fe', 'H'))[0]
        self.alpha_Fe = get_abnd_ratios(self.Zstar, ('alpha', 'Fe'))[0]
        self.Ca_Fe = get_abnd_ratios(self.Zstar, ('Ca', 'Fe'))[0]
        self.Si_Fe = get_abnd_ratios(self.Zstar, ('Si', 'Fe'))[0]
        self.Mg_Fe = get_abnd_ratios(self.Zstar, ('Mg', 'Fe'))[0]
        self.n=len(mstar) # number of star particles
        self.arrays = rec.array([idstar, mstar, self.pos[:,0], self.pos[:,1], self.pos[:,2], rstar, \
            self.vel[:,0], self.vel[:,1], self.vel[:,2], age, lookback_time, a, self.Zstar[:,0], self.Zstar[:,1], self.Zstar[:,2], \
            self.Zstar[:,3], self.Zstar[:,4], self.Zstar[:,5], self.Zstar[:,6], self.Zstar[:,7], self.Zstar[:,8], self.Zstar[:,9], \
            self.Zstar[:,10], np.log10(self.Zstar[:,0] / self.Zsol)],
            dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'),  ('r', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8'), ('age', '<f8'), \
            ('lookback_time','<f8'), ('a', '<f8'),('Z', '<f8'), ('He', '<f8'), ('C', '<f8'), ('N', '<f8'), ('O', '<f8'), ('Ne', '<f8'), \
            ('Mg', '<f8'), ('Si', '<f8'), ('S', '<f8'), ('Ca', '<f8'), ('Fe', '<f8'), \
            ('Z_Zsol', '<f8')])
        self.booly=np.zeros(len(self.arrays.id)) == 0 # initialize to all true, used to restrict particles

            
    # Put all arrays in order of increasing distance from center (argsort r)    
    def orderp(self):
        order = np.argsort(self.arrays.r)
        self.arrays.r = self.arrays.r[order]
        self.arrays.id = self.arrays.id[order]
        self.arrays.mass = self.arrays.mass[order]
        self.arrays.x = self.arrays.x[order]
        self.arrays.y = self.arrays.y[order]
        self.arrays.z = self.arrays.z[order]
        self.arrays.velx = self.arrays.velx[order]
        self.arrays.vely = self.arrays.vely[order]
        self.arrays.velz = self.arrays.velz[order]
        self.arrays.age = self.arrays.age[order]
        self.arrays.lookback_time = self.arrays.lookback_time[order]
        self.arrays.Z = self.arrays.Z[order]
        self.arrays.He = self.arrays.He[order]
        self.arrays.C = self.arrays.C[order]
        self.arrays.N = self.arrays.N[order]
        self.arrays.O = self.arrays.O[order]
        self.arrays.Ne = self.arrays.Ne[order]
        self.arrays.Mg = self.arrays.Mg[order]
        self.arrays.Si = self.arrays.Si[order]
        self.arrays.S = self.arrays.S[order]
        self.arrays.Ca = self.arrays.Ca[order]
        self.arrays.Fe = self.arrays.Fe[order]
        self.booly = self.booly[order]
        

    def restrict_age(self, age_lo, age_hi):
        """
        This function restricts ONLY the gas particles to lie within the specified temp range
        The function reset_to_orig can also be used to remove this condition
        
        """
        self.booly = ((self.arrays.age > age_lo) & (self.arrays.age < age_hi))
        
        
    def get_clr(self):
        return (0,0,1)
        
    # Label for Stars
    def get_labels(self):
        return 'Stars'
        
            
    def plot_alpha_Fe_vs_Fe_H_hist(self, alph = 'Ca', binx = 200, biny = 200, lims = [-1,-1]):
        """
        This function plots a single [Fe/H] vs [alpha/Fe] 2d histogram
        
        
        """
        if alph == 'Ca':
            if lims[0] == -1: # case for no colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.Ca_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
            else: # here we want colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.Ca_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(vmin=lims[0], vmax=lims[1]),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
                plt.xlabel(r'$\rm [Fe/H]$', fontsize = 22)
                plt.ylabel(r'$\rm [Ca/Fe]$', fontsize = 22)
        else: #plot all alpha
            if lims[0] == -1: # case for no colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.alpha_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
            else: # here we want colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.alpha_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(vmin=lims[0], vmax=lims[1]),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
                plt.xlabel(r'$\rm [Fe/H]$', fontsize = 22)
                plt.ylabel(r'$\rm [alpha/Fe]$', fontsize = 22)
                
		# Plot average values
        avg_Fe_H = np.mean(self.Fe_H[self.booly])
        if alph == 'Ca':
            avg_alph_Fe = np.mean(self.Ca_Fe[self.booly])
        else:
            avg_alph_Fe = np.mean(self.alpha_Fe[self.booly])          
            
        print("<[Fe_H]>, <[alph/Fe]>")
        print([avg_Fe_H,avg_alph_Fe])
        plt.axvline(x = avg_Fe_H, ls = 'dashed', color = 'k')
        plt.axhline(y = avg_alph_Fe, ls = 'dotted', color = 'k')
        plt.title(r"$\rm z = %0.4f$"%(self.header['redshift']))
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(r'$\rm Mass$ $\rm (M_{\odot})$', fontsize = 22)




    def age_scatter(self, agemin = 0, agemax = 13.8,  D = 5, Fs = 8, cmap = 'jet'):
        """
        
        """
    
        age = 13.8 - self.arrays.age
        self.restrict_parts(D)
        x = self.arrays.x[self.booly][((age[self.booly] >= agemin)&(age[self.booly] <= agemax))]
        y = self.arrays.y[self.booly][((age[self.booly] >= agemin)&(age[self.booly] <= agemax))]
        if len(x) > 0:
            plt.figure(1, figsize=(Fs, Fs))
            plt.scatter(x, y, c = age[self.booly][((age[self.booly] >= agemin)&(age[self.booly] <= agemax))], marker = '.', edgecolors = 'None', cmap = cmap)	
            plt.xlim(-5,5)
            plt.ylim(-5,5)	
            cbar = plt.colorbar()
            plt.clim(0,14)
            cbar.ax.set_ylabel(r'$\rm Age \, (Gyr)$', fontsize = 22)
            plt.xlabel(r'$\rm X \, (kpc)$', fontsize = 22)
            plt.ylabel(r'$\rm Y \, (kpc)$', fontsize = 22)
        self.reset_to_orig(recenter = True)
        self.recenter()

    def age_grad(self, nbins = 10, eqninbin = True, lbins = True):
        """
        
        """
        
        self.restrict_parts(self.rvir/10.)
        age = self.arrays.lookback_time[self.booly]	
        r = self.arrays.r[self.booly]
        
        binavg_or_sum, midbins, varbins, rbins = bin_data(rad = r, tobin = age, nbins = nbins, eqnuminbin = True, logbins = lbins)
        
        plt.figure(1, figsize=(Fs, Fs))
        plt.plot(np.log10(midbins), binavg_or_sum, marker = '^', color = self.color, label = self.label)	
        plt.xlim(-1.5,0.7)
        plt.ylim(0,13.8)	
        plt.xlabel(r'$\rm log \, R \, (kpc)$', fontsize = 22)
        plt.ylabel(r'$\rm Age \, (Gyr)$', fontsize = 22)
        
    def Fe_grad(self, nbins = 10, eqninbin = True, lbins = True):
        """
        
        """
        
        self.restrict_parts(D)
        Fe_H = self.Fe_H[self.booly]
        r = self.arrays.r[self.booly]
        
        binavg_or_sum, midbins, varbins, rbins = bin_data(rad = r, tobin = Fe_H, nbins = nbins, eqninbin = True, logbins = lbins)
        
        plt.figure(1, figsize=(Fs, Fs))
        plt.plot(np.log10(midbins), binavg_or_sum, marker = '^', color = self.color, label = self.label)	
        plt.xlim(-1.5,0.7)
        plt.ylim(-4,-2)	
        plt.xlabel(r'$\rm log \, R \, (kpc)$', fontsize = 22)
        plt.ylabel(r'$\rm [Fe/H]$', fontsize = 22)


    def plot_alpha_2Dhist(self, r_ext = 50, h_ext = 50, bins = 500, cmap = 'bone', cmin = 0, topdown = True, lims = [-1,-1]):
        """
        This function plots a visualization of the particle type from the z axis looking down,
            weighted by the alpha elements
        
        
        """
#         fig = plt.figure(1, figsize = (10,10))
#         ax1 = fig.add_subplot(1,1,1)
#         ax1.patch.set_facecolor('black')
        
        if topdown == True:
        	rr = np.sqrt((self.arrays.x[self.booly] **2) + (self.arrays.y[self.booly] **2))
        else:
        	rr = np.sqrt((self.arrays.x[self.booly] **2) + (self.arrays.z[self.booly] **2))

        
        # Only choose those particles in selected region
        x = self.arrays.x[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]
        y = self.arrays.y[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]
        z = self.arrays.z[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]

        weights = self.Ca_Fe[self.booly][((rr < r_ext) & \
            (np.abs(self.arrays.z[self.booly]) < h_ext))]
        
#         plot_m_2Dhist(x,y,z,weights,bins,cmap,cmin,topdown)
        if lims[0] == -1: # case for no colorbar scaling
        	plt.hist2d(x, y, bins = bins, norm=matplotlib.colors.LogNorm(), weights = weights)
        else: # here we want colorbar scaling
            plt.hist2d(x, y, bins = bins, norm=matplotlib.colors.LogNorm(vmin=lims[0], vmax=lims[1]), weights = weights)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(r'$\rm [Ca/Fe]$', fontsize = 22)

    def get_Mv_Lv(self, MvSun = 4.85, m2l = 1.0, r = 5.0):
        """
	
        Function for absolute magnitude of galaxy
        MvSun: absolute mag of sun
        m2l: stellar mass to light ratio
        r: radius of galaxy
	
        """
        Lv = m2l* self.calc_mass_within_r(rad = r, recenter = True)
        Mv = 4.85 - 2.5*np.log10(Lv)
        return Mv,Lv
    
    def plot_sfh(self, norm = True, gt = False, lw = 2, ls = 'solid', axes = True, Fs = 14, fntsz = 16, axwidth = 3, axlength = 12, color = None, label = None):
        """
        Takes an axis and plots the cumulative star formation history
        Tries to do this all internal, no need for other fcts
        
        """
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(Om0=self.header['omega_matter'], H0=self.header['hubble']*100)
        
        if color is not None:
        	clr = color
        else:
        	clr = self.color

        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1,1,1)
        ax1 = plt.gca()
        axwidth = axwidth
    	axlength = axlength
    	fntsz=fntsz
    	plt.rc('axes',linewidth=axwidth)
#     
#         self.recenter()
#         
#         self.restrict_parts(self.rvir/10.)
        arr = self.arrays.lookback_time[self.booly]		
        mpart = self.arrays.mass[self.booly]	
        ylog =  not norm

    	N = len(arr)
    	Xax = arr[np.argsort(arr)]
    	if norm == True: # norm the distrib
            Yax = [(N - i) / float(N) for i in range(N)]
    	else: # not normed
            masssort = mpart[np.argsort(arr)]
            Yax = np.cumsum(masssort)[::-1]
        if ylog == True:
            Yax = np.insert(Yax,0,Yax[0])
            Xax = np.insert(Xax,0,0)
            if label is None:
            	plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr, label = self.label)
            else:
            	plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr, label = label)
        else:
        	if label is None:
                    plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr, label = self.label)
                else:
                    plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr, label = label)

        plt.xlim(13.8,0)
        if ylog == False:
        	plt.ylim(0,1)
        else:
        	plt.ylim(2,6.5)
        
        if axes == True:
        
			zs = np.array([10,5,3,2,1,0.5,0.25,0])
			aa = [cosmo.lookback_time(z).value for z in zs]
	
			ax2 = ax1.twiny()
			ax2.set_xticks(aa)
			ax2.set_xticklabels(['{:g}'.format(z) for z in zs])
			ax2.set_xlim(ax1.get_xlim())
			zstr = ['$'+str(z)+'$' for z in zs]
			ax2.set_xticklabels(zstr, fontsize=fntsz)
			ax2.set_xlabel(r'$\rm Redshift$', fontsize = fntsz)
		
		
			xx = [12,10,8,6,4,2,0]
			xstr = ['$'+str(x)+'$' for x in xx]
			ax1.set_xticklabels(xstr[::-1], fontsize = fntsz)
			if norm == True:
				yy = [0.0,0.2,0.4,0.6,0.8,1.0]
				ystr = ['$'+str(y)+'$' for y in yy]
			else:
				yy = [2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5]
				ystr = ['$10^{'+str(y)+'}$' for y in yy]
			
			ax1.set_yticklabels(ystr, fontsize = fntsz)
			plt.xticks(weight = 'bold')
			plt.yticks(weight = 'bold')
		
			if norm == True:		
				ax1.set_ylabel(r'$\rm Cumulative \, Fractional \, SFH$', fontsize = fntsz)
			else:
				ax1.set_ylabel(r'$\rm M_{\star} \, (M_{\odot})$', fontsize = fntsz)
			ax1.set_xlabel(r'$\rm Lookback \, Time \, (Gyr)$', fontsize = fntsz, labelpad = 10)
				
			ax = ax1
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for line in ax.get_xticklines() + ax.get_yticklines():
				line.set_markersize(18)
				line.set_markeredgewidth(3)
			for tick in ax.xaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)
			for tick in ax.yaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)

			ax.tick_params(which='major',width=axwidth,length=axlength+5, pad = 10)
			ax.tick_params(which='minor',width=axwidth,length=axlength, pad  = 10)

			#####################

			ax = ax2
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for line in ax.get_xticklines() + ax.get_yticklines():
				line.set_markersize(18)
				line.set_markeredgewidth(3)
			for tick in ax.xaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)
			for tick in ax.yaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)

			ax.tick_params(which='major',width=axwidth,length=axlength+5, pad = 10)
			ax.tick_params(which='minor',width=axwidth,length=axlength, pad = 10)
			
			
			
			
			
			
    def plot_sfh_TEST(self, norm = True, gt = False, lw = 2, ls = 'solid', axes = True, Fs = 14, fntsz = 18, color = None, label = None):
        """
        Takes an axis and plots the cumulative star formation history
        Tries to do this all internal, no need for other fcts
        
        """
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(Om0=self.header['omega_matter'], H0=self.header['hubble']*100)
        
        if color is not None:
        	clr = color
        else:
        	clr = self.color

        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(1,1,1)
        ax1 = plt.gca()
        axwidth = 3
    	axlength = 12
    	fntsz=fntsz
    	plt.rc('axes',linewidth=axwidth)
#     
#         self.recenter()
#         
#         self.restrict_parts(self.rvir/10.)
        arr = self.arrays.lookback_time[self.booly]		
        mpart = self.arrays.mass[self.booly]	
        ylog =  not norm

    	N = len(arr)
    	Xax = arr[np.argsort(arr)]
    	if norm == True: # norm the distrib
            Yax = [(N - i) / float(N) for i in range(N)]
    	else: # not normed
            masssort = mpart[np.argsort(arr)]
            Yax = np.cumsum(masssort)[::-1]
        if ylog == True:
            Yax = np.insert(Yax,0,Yax[0])
            Xax = np.insert(Xax,0,0)
            if label is None:
            	plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr, label = self.label)
            else:
            	plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr, label = label)
        else:
        	if label is None:
                    plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr, label = self.label)
                else:
                    plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr, label = label)

        plt.xlim(13.8,12)
        if ylog == False:
        	plt.ylim(0,1)
        else:
        	plt.ylim(2,6.5)
        
        if axes == True:
        
			zs = np.array([25,20,15,10,5])
			aa = [cosmo.lookback_time(z).value for z in zs]
	
			ax2 = ax1.twiny()
			ax2.set_xticks(aa)
			ax2.set_xticklabels(['{:g}'.format(z) for z in zs])
			ax2.set_xlim(ax1.get_xlim())
			zstr = ['$'+str(z)+'$' for z in zs]
			ax2.set_xticklabels(zstr, fontsize=fntsz)
			ax2.set_xlabel(r'$\rm z$', fontsize = fntsz)
		
		
			xx = [13.5,13,12.5,12]
			xstr = ['$'+str(x)+'$' for x in xx]
			ax1.set_xticklabels(xstr[::-1], fontsize = fntsz)
			if norm == True:
				yy = [0.0,0.2,0.4,0.6,0.8,1.0]
				ystr = ['$'+str(y)+'$' for y in yy]
			else:
				yy = [2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5]
				ystr = ['$10^{'+str(y)+'}$' for y in yy]
			
			ax1.set_yticklabels(ystr, fontsize = fntsz)
			plt.xticks(weight = 'bold')
			plt.yticks(weight = 'bold')
		
			if norm == True:		
				ax1.set_ylabel(r'$\rm Cumulative \, Fractional \, SFH$', fontsize = fntsz)
			else:
				ax1.set_ylabel(r'$\rm M_{\star} \, (M_{\odot})$', fontsize = fntsz)
			ax1.set_xlabel(r'$\rm Lookback \, Time \, (Gyr)$', fontsize = fntsz)
				
			ax = ax1
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for line in ax.get_xticklines() + ax.get_yticklines():
				line.set_markersize(18)
				line.set_markeredgewidth(3)
			for tick in ax.xaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)
			for tick in ax.yaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)

			ax.tick_params(which='major',width=axwidth,length=axlength+5, pad = 10)
			ax.tick_params(which='minor',width=axwidth,length=axlength, pad = 10)

			#####################

			ax = ax2
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(fntsz)
			for line in ax.get_xticklines() + ax.get_yticklines():
				line.set_markersize(18)
				line.set_markeredgewidth(3)
			for tick in ax.xaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)
			for tick in ax.yaxis.get_minor_ticks():
				tick.label1.set_fontsize(fntsz/2)

			ax.tick_params(which='major',width=axwidth,length=axlength+5, pad = 10)
			ax.tick_params(which='minor',width=axwidth,length=axlength, pqd = 10)

    def plot_sfr(self, smooth, clr = None, specific = True, lw = 1, ls = '-', axes = True, \
    	Fs = 8, axwidth = 3, axlength = 12, fntsz = 22):
        """
        Plots the SFR smoothed over timescale smooth
        
        """
        
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(Om0=self.header['omega_matter'], H0=self.header['hubble']*100)
        
        fig = plt.figure(1, figsize=(Fs,Fs*0.75))
        ax1 = fig.add_subplot(1,1,1)
        ax1 = plt.gca()
    	plt.rc('axes',linewidth=axwidth)
        
        if clr is None:
            clr = self.color

        ttot = 13.8
        nbins = np.int(np.floor(ttot / (smooth/1000.)))
        bins = np.linspace(0,13.8,nbins)
        binavg_or_sum, midbins, varbins, rbins = bin_data(self.arrays.age[self.booly], \
                                                              self.arrays.mass[self.booly], nbins, logbins = False, suminbin = True)
        peryr = binavg_or_sum / (ttot*1e9/nbins)
        mtot = np.ndarray(nbins)
        mtot[0] = 0# np.sum(disk.arrays.mass)
        for i in range(nbins - 1):
            t = i+1
            mtot[t] = mtot[t-1] + binavg_or_sum[t]
        if smooth < 1000:
        	label = r"$\rm %i \, Myr \, Avg$"%(smooth)
        else:
        	label = r"$\rm %i \, Gyr \, Avg$"%(smooth/1000)
        if specific == True:
            logperyr = np.log10(peryr/mtot)
        else:
            logperyr = np.log10(peryr)
        logperyr[logperyr == -np.inf] = -10
        if specific == True:
            plt.plot(midbins, logperyr[::-1], label = label, color = clr, lw = lw, ls = ls)
            plt.ylabel(r"$\rm log \,sSFR \, [yr^{-1}]$", fontsize = fntsz,labelpad = 10)
        else:
            plt.plot(midbins, logperyr[::-1], label = label, color = clr, lw = lw, ls = ls)
            plt.ylabel(r"$\rm log \, SFR \, [M_* \, / \, yr^{-1}]$", fontsize = fntsz,labelpad = 10)
            
        plt.legend(loc = 1, numpoints = 1, frameon = False, fontsize = fntsz)
        plt.ylim(-5, -2)
        plt.xlim(13.8,0)
        
        if axes == True:
        
            zs = np.array([10,5,3,2,1,0.5,0.25,0])
            aa = [cosmo.lookback_time(z).value for z in zs]
            
            ax2 = ax1.twiny()
            ax2.set_xticks(aa)
            ax2.set_xticklabels(['{:g}'.format(z) for z in zs])
            ax2.set_xlim(ax1.get_xlim())
            zstr = ['$'+str(z)+'$' for z in zs]
            ax2.set_xticklabels(zstr, fontsize=fntsz*0.9)
            ax2.set_xlabel(r'$\rm Redshift$', fontsize = fntsz*1.2)
		
            xx = [12,10,8,6,4,2,0]
            xstr = ['$'+str(x)+'$' for x in xx]
            ax1.set_xticklabels(xstr[::-1], fontsize = fntsz)
            
            yy = [-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0]
            ystr = ['$'+str(y)+'$' for y in yy]
            
            ax1.set_yticklabels(ystr, fontsize = fntsz)
            plt.xticks(weight = 'bold')
            plt.yticks(weight = 'bold')
            
            
            ax1.set_ylabel(r'$\rm log \, SFR \, (M_{\odot} yr^{-1})$', fontsize = fntsz*1.2,labelpad = 10)
            ax1.set_xlabel(r"$\rm Lookback \, Time \, (Gyr)$", fontsize = fntsz*1.2, labelpad = 10)
            
            ax = ax1
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fntsz)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fntsz)
            for line in ax.get_xticklines() + ax.get_yticklines():
                line.set_markersize(18)
                line.set_markeredgewidth(3)
            for tick in ax.xaxis.get_minor_ticks():
                tick.label1.set_fontsize(fntsz/2)
            for tick in ax.yaxis.get_minor_ticks():
                tick.label1.set_fontsize(fntsz/2)
                
            ax.tick_params(which='major',width=axwidth,length=axlength+5, pad = 10)
            ax.tick_params(which='minor',width=axwidth,length=axlength, pad = 10)

			#####################

            ax = ax2
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(fntsz)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(fntsz)
            for line in ax.get_xticklines() + ax.get_yticklines():
                line.set_markersize(18)
                line.set_markeredgewidth(3)
            for tick in ax.xaxis.get_minor_ticks():
                tick.label1.set_fontsize(fntsz/2)
            for tick in ax.yaxis.get_minor_ticks():
                tick.label1.set_fontsize(fntsz/2)
                
            ax.tick_params(which='major',width=axwidth,length=axlength+5, pad = 10)
            ax.tick_params(which='minor',width=axwidth,length=axlength, pad = 10)


        return peryr








#############################################
class GasParticles(Particles):
    def get_arrays_from_file(self):
        rho=self.part['gas']['density'][:]
        idgas=self.part['gas']['id'][:]
        temp=self.part['gas']['temperature'][:]
        mgas=self.part['gas']['mass'][:]
        ne=self.part['gas']['electron.fraction'][:]
        nh=self.part['gas']['hydrogen.neutral.fraction'][:]
        self.pos=self.part['gas']['position'][:]
        sfr=self.part['gas']['sfr'][:]
        hsml=self.part['gas']['smooth.length'][:]
        self.vel=self.part['gas']['velocity'][:]
        rgas = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.Zgas=self.part['gas']['massfraction'][:]
        self.Fe_H = get_abnd_ratios(self.Zgas, ('Fe', 'H'))[0]
        self.alpha_Fe = get_abnd_ratios(self.Zgas, ('alpha', 'Fe'))[0]
        self.Ca_Fe = get_abnd_ratios(self.Zgas, ('Ca', 'Fe'))[0]
        self.n=len(mgas) # number of gas particles
        self.arrays = rec.array([idgas, mgas, self.pos[:,0], self.pos[:,1], self.pos[:,2], rgas, \
            self.vel[:,0], self.vel[:,1], self.vel[:,2], temp, rho, ne, nh, sfr, hsml, self.Zgas[:,0], \
            self.Zgas[:,1], self.Zgas[:,2], self.Zgas[:,3], self.Zgas[:,4], self.Zgas[:,5], self.Zgas[:,6], self.Zgas[:,7], \
            self.Zgas[:,8], self.Zgas[:,9], self.Zgas[:,10], np.log10(self.Zgas[:,0] / self.Zsol)], \
            dtype = [('id', '<i8'), ('mass', '<f8'), \
            ('x', '<f8'), ('y', '<f8'), ('z', '<f8'), ('r', '<f8'),  ('velx', '<f8'), ('vely', '<f8'), \
            ('velz', '<f8'), ('temp', '<f8'), ('rho', '<f8'), ('ne', '<f8'), \
            ('nh', '<f8'), ('sfr', '<f8'), ('hsml', '<f8'), ('Z', '<f8'), ('He', '<f8'), \
            ('C', '<f8'), ('N', '<f8'), ('O', '<f8'), ('Ne', '<f8'), ('Mg', '<f8'), \
            ('Si', '<f8'), ('S', '<f8'), ('Ca', '<f8'), ('Fe', '<f8'), \
            ('Z_Zsol', '<f8')])
        self.booly=np.zeros(len(self.arrays.id)) == 0 # initialize to all true, used to restrict particles


    # Put all arrays in order of increasing distance from center (argsort r)    
    def orderp(self):
        order = np.argsort(self.arrays.r)
        self.arrays.r = self.arrays.r[order]
        self.arrays.id = self.arrays.id[order]
        self.arrays.mass = self.arrays.mass[order]
        self.arrays.x = self.arrays.x[order]
        self.arrays.y = self.arrays.y[order]
        self.arrays.z = self.arrays.z[order]
        self.arrays.velx = self.arrays.velx[order]
        self.arrays.vely = self.arrays.vely[order]
        self.arrays.velz = self.arrays.velz[order]
        self.arrays.temp = self.arrays.temp[order]
        self.arrays.rho = self.arrays.rho[order]
        self.arrays.ne = self.arrays.ne[order]
        self.arrays.nh = self.arrays.nh[order]
        self.arrays.sfr = self.arrays.sfr[order]
        self.arrays.hsml = self.arrays.hsml[order]
        self.arrays.Z = self.arrays.Z[order]
        self.arrays.He = self.arrays.He[order]
        self.arrays.C = self.arrays.C[order]
        self.arrays.N = self.arrays.N[order]
        self.arrays.O = self.arrays.O[order]
        self.arrays.Ne = self.arrays.Ne[order]
        self.arrays.Mg = self.arrays.Mg[order]
        self.arrays.Si = self.arrays.Si[order]
        self.arrays.S = self.arrays.S[order]
        self.arrays.Ca = self.arrays.Ca[order]
        self.arrays.Fe = self.arrays.Fe[order]
        self.booly = self.booly[order]

    # Color for gas is always red
    def get_clr(self):
        return (1,0,0)
        
    # Label for Gas
    def get_labels(self):
        return 'Gas'
    
    def restrict_temp(self, temp_lo, temp_hi):
        """
        This function restricts ONLY the gas particles to lie within the specified temp range
        The function reset_to_orig can also be used to remove this condition
        
        """
        self.booly[((self.arrays.temp < temp_lo) | (self.arrays.temp > temp_hi))] = False
            
    def plot_alpha_Fe_vs_Fe_H_hist(self, alph = 'Ca', binx = 200, biny = 200, lims = [-1,-1]):
        """
        This function plots a single [Fe/H] vs [alpha/Fe] 2d histogram
        
        
        """
        if alph == 'Ca':
            if lims[0] == -1: # case for no colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.Ca_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
            else: # here we want colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.Ca_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(vmin=lims[0], vmax=lims[1]),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
                plt.xlabel(r'$\rm [Fe/H]$', fontsize = 22)
                plt.ylabel(r'$\rm [Ca/Fe]$', fontsize = 22)
        else: #plot all alpha
            if lims[0] == -1: # case for no colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.alpha_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
            else: # here we want colorbar scaling
                plt.hist2d(self.Fe_H[self.booly], self.alpha_Fe[self.booly], \
                               norm=matplotlib.colors.LogNorm(vmin=lims[0], vmax=lims[1]),\
                               bins = [binx, biny], weights = self.arrays.mass[self.booly])
                plt.xlabel(r'$\rm [Fe/H]$', fontsize = 22)
                plt.ylabel(r'$\rm [alpha/Fe]$', fontsize = 22)
                
		# Plot average values
        avg_Fe_H = np.mean(self.Fe_H[self.booly])
        if alph == 'Ca':
            avg_alph_Fe = np.mean(self.Ca_Fe[self.booly])
        else:
            avg_alph_Fe = np.mean(self.alpha_Fe[self.booly])
            
        print("<[Fe_H]>, <[alph/Fe]>")
        print([avg_Fe_H,avg_alph_Fe])
        plt.axvline(x = avg_Fe_H, ls = 'dashed', color = 'k')
        plt.axhline(y = avg_alph_Fe, ls = 'dotted', color = 'k')
        plt.title(r"$\rm z = %0.4f$"%(self.header['redshift']))
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(r'$\rm Mass$ $\rm (M_{\odot})$', fontsize = 22)


    def plot_T_vs_density(self, clr = 'k'):
        """
        
        """
        plt.figure(1, figsize=(10,10))
        plt.plot(np.log10(self.arrays.rho[self.booly]), np.log10(self.arrays.temp[self.booly]), ls = ' ', marker = '.', color = clr) 
        plt.ylabel(r'$\rm log \, Temperature \, (K)$', fontsize = 22)
        plt.xlabel(r'$\rm log \, Density \, (M_{sun} / kpc^3)$', fontsize = 22)
        plt.title(r"$\rm z = %0.4f$"%(self.header['redshift']))

    def plot_T_vs_density_hist(self, binx = 200, biny = 200, lims = [-1,-1]):
        """
        
        
        """
        plt.figure(1, figsize=(10,10))
        if lims[0] == -1: # case for no colorbar scaling
            plt.hist2d(np.log10(self.arrays.rho[self.booly]), np.log10(self.arrays.temp[self.booly]), \
                norm=matplotlib.colors.LogNorm(),\
                bins = [binx, biny], weights = self.arrays.mass[self.booly])
        else: # here we want colorbar scaling
            plt.hist2d(np.log10(self.arrays.rho[self.booly]), np.log10(self.arrays.temp[self.booly]), \
                norm=matplotlib.colors.LogNorm(vmin=lims[0], vmax=lims[1]),\
                bins = [binx, biny], weights = self.arrays.mass[self.booly])
        plt.ylabel(r'$\rm log \, Temperature \, (K)$', fontsize = 22)
        plt.xlabel(r'$\rm log \, Density \, (cm^3)$', fontsize = 22)
        plt.title(r"$\rm z = %0.4f$"%(self.header['redshift']))
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(r'$\rm Mass$ $\rm (M_{\odot})$', fontsize = 22)
#         plt.xlim(3,8)
#         plt.ylim(2,7)
        

 
   
class DiskParticles(Particles):
    def get_arrays_from_file(self):
        idstar=self.part['disk']['id'][:] # the ids of the star particles
        mstar=self.part['disk']['mass'][:] * 10 ** 10 / self.hh
        self.pos=self.part['disk']['position'][:]/ self.hh
        self.vel=self.part['disk']['velocity'][:]
        rstar = np.sqrt(self.pos[:,0]**2+self.pos[:,1]**2+self.pos[:,2]**2)
        self.n=len(mstar) # number of disk particles
        self.arrays = rec.array([idstar, mstar, self.pos[:,0], self.pos[:,1], self.pos[:,2], rstar, \
            self.vel[:,0], self.vel[:,1], self.vel[:,2]],
            dtype = [('id', '<i8'), ('mass', '<f8'), ('x', '<f8'), ('y', '<f8'), \
            ('z', '<f8'),  ('r', '<f8'), ('velx', '<f8'), ('vely', '<f8'), ('velz', '<f8')])
        self.booly=np.zeros(len(self.arrays.id)) == 0 # initialize to all true, used to restrict particles
        

# ########################################################################################
# Scrips that aren't methods

def bin_data(rad = R_proj, tobin = vz, nbins = 11, logbins = False, eqnuminbin = False, suminbin = False, allspace = True):
    """
    Takes an array and a number of bins and returns the bin midpoints.
    
    """
    if logbins == True:
        rad = np.log10(rad)
        
    # Space equally (linear or log)
    if eqnuminbin == False:
    	if allspace == False:
        	rmax = np.amax(rad) # largest radius
        else:
        	rmax = 13.8
        rmin = np.amin(rad) # smallest radius
        rbins = np.linspace(rmin, rmax, nbins+1) # bins in r
    
    # Equal number in each bin
    else:
        rbins = np.ndarray(nbins+1)
        for i in range(nbins+1):
            rbins[i] = stats.scoreatpercentile(rad, i * 100. / nbins)    
    
    if logbins == True:
        rbins = 10 ** rbins
        
    midbins = np.ndarray(nbins)
    varbins = np.ndarray(nbins)
    binavg_or_sum = np.ndarray(nbins)
    ninbin = np.zeros(nbins)
    for i in range(nbins):
        midbins[i] = (rbins[i] + rbins[i+1]) / 2.
        # Now split data into bins and return average in bin or sum in bin
        stuff_in_bin = tobin[((rad > rbins[i]) & (rad < rbins[i+1]))]
        ninbin[i] = len(stuff_in_bin)
        if suminbin == True:
            binavg_or_sum[i] = np.sum(stuff_in_bin)
        else: # want average in bin
            binavg_or_sum[i] = np.sum(stuff_in_bin) / len(stuff_in_bin)
        varbins[i] = np.var(stuff_in_bin)
    
    return binavg_or_sum, midbins, varbins, rbins
    
    
def get_all(model, runnum, snapstart, snapend, snapskip):
    """
    Used for a list of snapshots, given by first, last and number between each,
    to get arrays of snaps of each particle type
    
    """
    snaplist = np.arange(snapstart, snapend+1, snapskip) # could be single or list, always array
    snapfulls = []
    snaps = []
    fnames = []
    arestars = []
    darkp = []
    gasp = []
    starsp = []


    for snapdig in snaplist:
        # Get full snap list
        snap = '{:04d}'.format(snapdig)
        snaps.append(snap) # Just the threee digits
        fname = "discus_igni_easy_%s.hdf5"%(snap)
        fnames.append(fname)
        if header['particle.numbers.in.file'][0] == 0: # no stars case
            isstars = bool(False)
        else:
            isstars = bool(True)
        arestars.append(isstars)
        dm = DarkParticles(fname, snapdig)
        gas = GasParticles(fname, snapdig)
        if isstars:
            stars = StarParticles(fname, snapdig)
        
        # Recenter based on location of ALL particles if star particle exists
        if isstars:
            pos = np.vstack((dm.pos, stars.pos, gas.pos))
            m = np.hstack((dm.arrays.mass, stars.arrays.mass, gas.arrays.mass))
            v = np.vstack((dm.vel, stars.vel, gas.vel))
        else: # Only use dark matter and gas if no stars
            pos = np.concatenate((dm.pos, gas.pos), axis=0)
            m = np.concatenate((dm.arrays.mass, gas.arrays.mass))
            v = np.concatenate((dm.vel, gas.vel), axis=0)

        dm.recenter(pos,v,m)
        gas.recenter(pos,v,m)
        if isstars:
            stars.recenter(pos,v,m)
        
        # Now append to list of snaps for each type of particle
        darkp.append(dm)
        gasp.append(gas)
        if isstars:
            starsp.append(stars)
        else:
            starsp.append(-1)
            
    return darkp, gasp, starsp, arestars

def order_all(dm, gas, stars):
    for t in range(len(dm)):
        dm[t].orderp()
        gas[t].orderp()
        if stars[t] != -1:
            stars[t].orderp() 

def density_profiles(dm, gas, nsnaps, model, runnum, xlims=[-3,2.3], ylims=[0,10]):
    """
    """
    for t in range(nsnaps):
        snap3 = '{0:03}'.format(dm[t].snap)
        outdir = "%s/%s/%s_%s/profiles"%(workpath, model, model, runnum)
        plt.figure(1, figsize=(10,10))
        dm[t].plot_rho()
        gas[t].plot_rho()
        outfile = "%s/dm_gas_density_snap_%s.jpg"%(outdir,snap3)
        plt.xlim(xlims[0],xlims[1])
        plt.ylim(ylims[0], ylims[1])
        #plt.title(r"$\rm %.2f \, Myr$"%(dm[t].snap*0.05), fontsize = 30)
        plt.title(r"$\rm %d0 \, Myr$"%(dm[t].snap), fontsize = 30)
        plt.savefig(outfile)
        plt.clf()
        
def get_abnd_ratios(Z, *args):

	""" 
	Return the specified metallicity from the simulation
	
	Parameters
	-----------
	
	Z: array of metallicites
	args: arguments should be either comma-separated tuples or a list of
			tuples, e.g. (spec1, spec2)
	     	spec1: species in numerator, str: H is not an option
			spec2: species in denominator, str: alpha is not option
	"""

	#Dictionaries listed above
	#spec_dict = dw_const.spec_dict
	#solar_abnds = dw_const.solar_abnds
	#amu = dw_const.amu
	
	spec_dict = {'Z': 0, 'H': 0, 'He': 1, 'C': 2, 'N': 3, 'O': 4, 'Ne': 5, \
		'Mg': 6, 'Si': 7, 'S': 8, 'Ca': 9, 'Fe': 10}
		
		#From Anders & Grevesse 1989, Fe from Sneden et al 1992
	solar_abnds = {'Mg': 7.58, 'Ca': 6.36, 'Si': 7.55, 'Ti': 4.99, 'Fe': 7.52, 'O': 8.93,\
			   'H': 12.}
			   
	#Atomic mass units
	amu = {'H': 1.008, 'He': 4.003, 'Mg': 24.305, 'Ca': 40.078, 'Si': 28.086, 'Ti': 47.88,\
	   'Fe': 55.933, 'O': 16.00}

	#Check of list, single, or multiple inputs, and make into list
	if len(args) == 1:
		if isinstance(args, list):
			args_list = args
		else:
			args_list = list(args)
	elif len(args) > 1:
		args_list = [a for a in args]
	else:
		sys.stderr.write('Must enter an argument')
		sys.exit()

	#Check if tuple
	for a in args_list:
		if isinstance(a, tuple):
			pass
		else:
			sys.stderr.write('Each argument must be a tuple')
			sys.exit()

	#Check if metallicity data exists for species
	elems = []
	for a in args_list:
		elems += [a[0]]
		elems += [a[1]]
	elems = list(set(elems))
	keys = spec_dict.keys()
	for elem in elems:
		if elem not in keys:
			if elem != 'alpha':
				sys.stderr.write('Invalid metallicity species {}\n'.format(elem))
				sys.exit()
	
	#Check if alpha/Fe in arguments
	if ('alpha', 'Fe') in args_list:
		if ('Mg', 'Fe') in args_list: pass
		else: args_list += [('Mg', 'Fe')]
		if ('Ca', 'Fe') in args_list: pass
		else: args_list += [('Ca', 'Fe')]
		if ('Si', 'Fe') in args_list: pass
		else: args_list += [('Si', 'Fe')]
		alpha = True
		#Remove from arguments
		index = [i for i,(t1,t2) in enumerate(args_list) if t1 == 'alpha']
		del args_list[index[0]]
	else:
		alpha = False

	#Get data to calculate abundance ratios
	
	abnd_ratios = []
	for a in args_list:
		spec1 = a[0]; spec2 = a[1]
		
	#Z is a 2-D array, where each element corresponds to the metallicites for an
	#individual star particle
	#Each star particle has 12 metallicity species

		#Get X1 and X2 metal mass fractions, solar abundances, and amus
		for spec in spec_dict.keys():

			if spec1 == spec:
				X1 = Z[:, spec_dict[spec1]]
				m_X1 = amu[spec1]
				solar_X1_abnd = solar_abnds[spec1]

			if spec2 == spec:
				X2_tmp = Z[:, spec_dict[spec2]]
				if spec2 == 'H':
					X2 = 1 - (X2_tmp + Z[:,1])
				else:
					X2 = X2_tmp 
				solar_X2_abnd = solar_abnds[spec2]
				m_X2 = amu[spec2]

		#Calculate abundance ratio
		solar_ratio = solar_X1_abnd - solar_X2_abnd
		abnd_ratio = np.log10(m_X2*X1/(m_X1*X2)) - solar_ratio
# 		abnd_ratio = m_X2*X1/(m_X1*X2) / (10**solar_ratio)
		abnd_ratios += [abnd_ratio]
		
	abnd_ratios = np.array(abnd_ratios)

	#Calculate average alpha/Fe ratio
	if alpha == True:
	
		ind = [i for i,t in enumerate(args_list) if (t == ('Mg', 'Fe') or t == ('Ca', 'Fe') or t == ('Si', 'Fe'))]
		ind = np.array(ind)
		alpha_abnd = np.sum(abnd_ratios[ind], axis=0)/float(len(abnd_ratios[ind]))
		abnd_ratios = np.insert(abnd_ratios, index[0], alpha_abnd, axis=0)
		
	return abnd_ratios


def plot_cumu_simple(arr, labels = 'XXX', clr = 'k', lw = 2, ls = 'solid', loc = 1, fntsz = 22, xlog=False, ylog=False, norm=True, gt=True):
    """
    Yes another fct to plot the cumulative distrib of some array, in this case N(< A)
    Here there are options for log-log plot, normed plot, greater than A or less than

    ARGUMENTS:
    arr - the array in consideration
    labels - optional label
    xlog - plot x in log space?
    ylog - plot y in log space?
    gt - True of plotting N(>arr)
    norm - True if normed distribution desired 
    
    """
    N = len(arr)
    Xax = arr[np.argsort(arr)]

    if gt == True: # plot N(>A)
    	if norm == True:
        	Yax = [(N - i) / float(N) for i in range(N)]
        else:
        	Yax =  [(N - i) for i in range(N)]
    else: # plot N(<A)
        if norm == True: 
        	Yax = [(i+1) / float(N) for i in range(N)]
        else:
        	Yax = [(i+1) for i in range(N)]

    if xlog == True:
        if ylog == True:
            if labels == 'XXX':
                plt.plot(np.log10(Xax), np.log10(Yax), lw = lw, ls = ls, color = clr)
            else:
                plt.plot(np.log10(Xax), np.log10(Yax), lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
        else:
            if labels == 'XXX':
                plt.plot(np.log10(Xax), Yax, lw = lw, ls = ls, color = clr)
            else:
                plt.plot(log10(Xax), Yax, lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
    else:
        if ylog == True:
            if labels == 'XXX':
                plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr)
            else:
                plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
        else:
            if labels == 'XXX':
                plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr)
            else:
                plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
#     plt.xlabel(r'$\rm Formation \, Time \, (Gyr)$', fontsize=22)    
    
#     if norm == True:
#     	plt.ylabel(r'$\rm frac(< \, Formation \, Time)$', fontsize=22)
#     else:
#     	plt.ylabel(r'$\rm log \, M_{*}(< \, Formation \, Time)$', fontsize=22)

def plot_cumu(arr, mass, labels = 'XXX', clr = 'k', lw = 2, ls = 'solid', fntsz = 22, xlog=False, ylog=False, norm=True, gt=True):
    """
    Yes another fct to plot the cumulative distrib of some array, in this case N(< A)
    Here there are options for log-log plot, normed plot, greater than A or less than

    ARGUMENTS:
    arr - the array in consideration
    mass -  array of masses
    labels - optional label
    xlog - plot x in log space?
    ylog - plot y in log space?
    gt - True of plotting N(>arr)
    norm - True if normed distribution desired 
    
    """
    N = len(arr)
    Xax = arr[np.argsort(arr)]
    if norm == True: # norm the distrib
        if gt == True: # plot N(>A)
            Yax = [(N - i) / float(N) for i in range(N)]
        else: # plot N(<A)
            Yax = [(i+1) / float(N) for i in range(N)]
    else: # not normed
    	masssort = mass[np.argsort(arr)]
        if gt == True: # plot M(>A)
        	Yax = np.cumsum(masssort)[::-1]
        else: # plot N(<A)
            Yax = np.cumsum(masssort)
    if xlog == True:
        if ylog == True:
            if labels == 'XXX':
                plt.plot(np.log10(Xax), np.log10(Yax), lw = lw, ls = ls, color = clr)
            else:
                plt.plot(np.log10(Xax), np.log10(Yax), lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
        else:
            if labels == 'XXX':
                plt.plot(np.log10(Xax), Yax, lw = lw, ls = ls, color = clr)
            else:
                plt.plot(log10(Xax), Yax, lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
    else:
        if ylog == True:
            if labels == 'XXX':
                plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr)
            else:
                plt.plot(Xax, np.log10(Yax), lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
        else:
            if labels == 'XXX':
                plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr)
            else:
                plt.plot(Xax, Yax, lw = lw, ls = ls, color = clr, label = labels)
                plt.legend(loc = loc, frameon = False, numpoints = 1, fontsize = fntsz)
#     plt.xlabel(r'$\rm Formation \, Time \, (Gyr)$', fontsize=22)    
	plt.xlabel(r'$\rm Temperature (K)$', fontsize=22)    

    plt.xlim(0,13.8)
    
    if norm == True:
    	plt.ylabel(r'$\rm frac(< \, Formation \, Time)$', fontsize=22)
    else:
    	plt.ylabel(r'$\rm log \, M_{*}(< \, Formation \, Time)$', fontsize=22)


def gen_bins(rad, nbins, logbins = False, eqnuminbin = False):
    """
    Takes an array and a number of bins and returns the bin edges.
    
    """
    if logbins == True:
        rad = np.log10(rad)
        
    # Space equally (linear or log)
    if eqnuminbin == False:
        rmax = np.amax(rad) # largest radius
        rmin = np.amin(rad) # smallest radius
        rbins = np.linspace(rmin, rmax, nbins+1) # bins in r
    
    # Equal number in each bin
    else:
        rbins = np.ndarray(nbins+1)
        for i in range(nbins+1):
            rbins[i] = stats.scoreatpercentile(rad, i * 100. / nbins)    
    
    if logbins == True:
        rbins = 10 ** rbins
        
    return rbins
    
def plot_m_2Dhist(x,y,z,weights,bins=1000,cmap = 'bone', cmin = 0, topdown = True):
	"""
	This function plots a visualization of the particle type from the z axis looking down,
		weighted by the particle mass
	
	
	"""
	fig = plt.figure(1, figsize = (10,10))
	ax1 = fig.add_subplot(1,1,1)
	ax1.patch.set_facecolor('black')
	if topdown == True:
		ax1.hist2d(x, y, norm=matplotlib.colors.LogNorm(), bins = bins, \
			weights = weights, cmap = cmap, cmin = cmin)
		plt.xlabel(r'$\rm X (kpc)$', fontsize = 22)
		plt.ylabel(r'$\rm Y (kpc)$', fontsize = 22)

	else: # side view
		ax1.hist2d(x, z, norm=matplotlib.colors.LogNorm(), bins = bins, \
			weights = weights, cmap = cmap, cmin = cmin)
		plt.xlabel(r'$\rm X (kpc)$', fontsize = 22)
		plt.ylabel(r'$\rm Z (kpc)$', fontsize = 22)


# def get_r_birth(stars, new_ids, snapstart, snapend):
#     """
#     This function takes a list of unique ids and searches all snaps for the first time it 
#         appears and gets the current (birth) radius
#     Returns snapnum, r_birth
#     """
# 
#     r_birth = np.array(len(new_ids))
#     #Now get r_birth
#     for id in new_ids:
#         nsnaps = 1+(snapend - snapstart) # Number of total snapshots
#         snaplist = np.arange(snapstart, snapend+1)
#         for snapdig in snaplist:
#             # Get full snap list
#             snap = '{:04d}'.format(snapdig)
#             fname = 'discus_igni_easy_%s.hdf5'%(snap)
#             part, header = readfire(snap)
#             dm = DarkParticles(fname, snapdig, part, header)
#             stars = StarParticles(fname, snapdig, part, header)
#             gas = GasParticles(fname, snapdig, part, header)
#             pos = np.vstack((dm.pos, stars.pos, gas.pos))
#             m = np.hstack((dm.arrays.mass, stars.arrays.mass, gas.arrays.mass))
#             v = np.vstack((dm.vel, stars.vel, gas.vel))
#         #     dm.recenter(pos,v,m)
#             stars.recenter(pos,v,m)
#         #     gas.recenter(pos,v,m)
#             if stars
#             r_birth[t] = stars.arrays.r(np.where(stars.arrays.id == new_ids[t])[0])
        
        
def compare_sfr(stars, smooth1 = 10, smooth2 = 200, fac = 1.5, specific = False, clr = 'k', lw = 1.5):
    """
    calculates fraction of mass formed in burst (smooth1 = fac * smooth2 avg)
    
    """
#     ttot = np.amax(stars.arrays.age) - np.amin(stars.arrays.age)
    ttot = 13.8
    nbins1 = np.int(np.floor(ttot / (smooth1/1000.)))
    nbins2 = np.int(np.floor(ttot / (smooth2/1000.)))
    binsum1, midbins1, varbins1, binedges1 = cf.bin_data(stars.arrays.age, stars.arrays.mass, nbins1, logbins = False, \
        suminbin = True)
    binsum2, midbins2, varbins2, binedges2 = cf.bin_data(stars.arrays.age, stars.arrays.mass, nbins2, logbins = False, \
        suminbin = True)
    peryr1 = binsum1 / (ttot*1e9/nbins1)
    peryr2 = binsum2 / (ttot*1e9/nbins2)
    
    mburst = np.ndarray(nbins2)
    mpostburst = np.ndarray(nbins2)
    msteady = np.ndarray(nbins2)
    for i in range(nbins2):
    	bin1_in_bin2 = peryr1[((midbins1 > binedges2[i]) & (midbins1 < binedges2[i+1]))]
    	mburst[i] = np.sum(bin1_in_bin2[bin1_in_bin2 > fac * peryr2[i]])
    	mpostburst[i] = np.sum(bin1_in_bin2[bin1_in_bin2 < peryr2[i]]/fac)
    	msteady[i] = np.sum(bin1_in_bin2[((bin1_in_bin2 < fac * peryr2[i]) & (bin1_in_bin2 > peryr2[i]/fac))])
		
    logmburst = np.log10(mburst)
    logmburst[logmburst == -np.inf] = -10
    logmpostburst = np.log10(mpostburst)
    logmpostburst[logmpostburst == -np.inf] = -10
    logmsteady = np.log10(msteady)
    logmsteady[logmsteady == -np.inf] = -10
    
    fburst = np.sum(mburst) / (np.sum(mpostburst)+np.sum(msteady)+np.sum(mburst))
    
    labelB = r"$\rm SFR_{%iMyr} > %.1f \times SFR_{%iMyr}$"%(smooth1,fac,smooth2)
    labelPB = r"$\rm SFR_{%iMyr} < SFR_{%iMyr} / %.1f$"%(smooth1,smooth2,fac)
    labelS = r"$\rm SFR_{%iMyr}/%.1f < SFR_{%iMyr} < %.1f \times SFR_{%iMyr}$"%(smooth2,fac,smooth1,fac,smooth2)
    plt.figure(1, figsize=(10,10))
    plt.plot(midbins2, logmburst, label = labelB, color = clr, lw = lw*1.5)
    plt.plot(midbins2, logmpostburst, label = labelPB, color = 'k', lw = lw, ls = 'dotted')
    plt.plot(midbins2, logmsteady, label = labelS, color = 'grey', lw = lw)
    plt.ylabel(r"$\rm (M_{\odot} \, yr^{-1})$", fontsize = 18)
    
    plt.xlabel(r"$\rm Formation \, Time \, (Gyr)$", fontsize = 18)
    plt.legend(loc = 1, numpoints = 1, frameon = False, fontsize = 16)
    plt.ylim(-5, -1)
    plt.xlim(0,13.8)

def compare_sfr_sig(stars, smooth1 = 10, smooth2 = 100):
    """
    calculates variance between sfr smoothed over smooth1 < smooth2 (units = Myr)
    
    """
#     ttot = np.amax(stars.arrays.age) - np.amin(stars.arrays.age)
    ttot = 13.8
    # number of bins = age universe / years in bins - 1000 accounts for Myr vs Gyr
    nbins1 = np.int(np.floor(ttot / (smooth1/1000.)))
    nbins2 = np.int(np.floor(ttot / (smooth2/1000.)))
    # bin up both time and mass - varbins not used
    binsum1, midbins1, varbins1, binedges1 = bin_data(stars.arrays.age, stars.arrays.mass, nbins1, logbins = False, suminbin = True)
    binsum2, midbins2, varbins2, binedges2 = bin_data(stars.arrays.age, stars.arrays.mass, nbins2, logbins = False, \
        suminbin = True)
        
    # calculate SFR in each bin, over each timescale (smooth1 and smooth2) - 1e9 accounts for Gyr (to get Msun/yr)
    peryr1 = binsum1 / (ttot*1e9/nbins1)
    peryr2 = binsum2 / (ttot*1e9/nbins2)
    
    # Make an array of length nbins1 where each entry is a duplicate of the value for that time in bins2
    ratbins = nbins1/nbins2
    peryr2_ext = np.repeat(peryr2,ratbins)
 
 	# Now calculate for each small bin, the ratio of the SFR formed in that bin to the average over the larger bin
    rat = peryr1/peryr2_ext
    return np.std(np.log(rat[rat > 1e-07])) # don't take dispersion of inf
 
 	
 
 
 
def age_scatter(stars, D, agemin = 0, agemax = 13.8, Fs = 8, cmap = 'jet'):
	"""
	
	"""
 	age = 13.8 - stars.arrays.age
	stars.restrict_parts(D)
	x = stars.arrays.x[stars.booly][((age[stars.booly] >= agemin)&(age[stars.booly] <= agemax))]
	y = stars.arrays.y[stars.booly][((age[stars.booly] >= agemin)&(age[stars.booly] <= agemax))]
	if len(x) > 0:
		plt.figure(1, figsize=(Fs, Fs))
		plt.scatter(x, y, c = age[stars.booly][((age[stars.booly] >= agemin)&(age[stars.booly] <= agemax))], marker = '.', edgecolors = 'None', cmap = cmap)	
		plt.xlim(-10,10)
		plt.ylim(-10,10)	
		cbar = plt.colorbar()
		plt.clim(0,14)
		cbar.ax.set_ylabel(r'$\rm Age \, (Gyr)$', fontsize = 22)
		plt.xlabel(r'$\rm X \, (kpc)$', fontsize = 22)
		plt.ylabel(r'$\rm Y \, (kpc)$', fontsize = 22)
 
 

def plot_Lv_rhalf(jdata, jdataM, desdata, des2data):
	"""
	
	"""
	import sat_fct as sf
	import all_fcts as afct
	Fs = 14
	
# 	jdata = open('dsph_ppties_labels.dat', 'r')
	jlines = jdata.readlines()
	# jdataM = open('dsph_mstar_labels.dat', 'r')
	jlinesM = jdataM.readlines()

# 	desdata = open('DES_dwfs_ppties.dat', 'r')
	deslines = desdata.readlines()

	# Year 2 data
# 	des2data = open('DES_Y2_dwfs_ppties.dat', 'r')
	des2lines = des2data.readlines()

	to3D = 1 # If using 3D, set to 4./3, 2D set to 1
	to2D = 3./4 # If using 3D, set to 1, 2D set to 3./4, 

	# ## Everything is in 3D, so we multiply Re_2d by 4/3
	axwidth = 3
	axlength = 12
	fntsz=32
	plt.rc('axes',linewidth=axwidth)
	fig = plt.figure(1, figsize=(Fs, Fs))
	ax = plt.gca()
	Name, RhalfDES, er_pl, er_mn, MstarDES, em_pl, em_mn = sf.get_desdata(deslines)
	RhalfDES = RhalfDES * to3D
	er_pl = er_pl * to3D
	er_mn = er_mn * to3D
	MstarDES = MstarDES * 1e3
	em_pl = em_pl * 1e3
	em_mn = em_mn * 1e3

	# Year 2 data
	Name2, RhalfDES2, er_pl2, er_mn2, MstarDES2, em_pl2, em_mn2 = sf.get_desdata(des2lines)
	RhalfDES2 = RhalfDES2 * to3D
	er_pl2 = er_pl2 * to3D
	er_mn2 = er_mn2 * to3D
	MstarDES2 = MstarDES2 * 1e3
	em_pl2 = em_pl2 * 1e3
	em_mn2 = em_mn2 * 1e3

	plt.plot(MstarDES, RhalfDES, marker = 'h', color = 'k', ls = ' ', markersize = 6)
	plt.errorbar(MstarDES, RhalfDES, xerr = [em_mn, em_pl], yerr = [er_mn, er_pl], color = 'k', ls = ' ')

	plt.plot(MstarDES2, RhalfDES2, marker = '.', color = 'k', ls = ' ', markersize = 18)
	plt.errorbar(MstarDES2, RhalfDES2, xerr = [em_mn2, em_pl2], yerr = [er_mn2, er_pl2], color = 'k', ls = ' ')

	#### Here is data for Hydra II (Martin et al. 2015)
	msHydr = 7.94 * 1e3
	em_pl_H = 2.1 * 1e3
	em_mn_H = 1.6 * 1e3
	RhHydr = 68
	er_pl_H = 11
	er_mn_H = 11

	plt.plot(msHydr, RhHydr, marker = 's', color = 'k', ls = ' ', markersize = 10)
	plt.errorbar(msHydr, RhHydr, xerr = np.array([[em_mn_H], [em_pl_H]]), yerr = er_mn_H, color = 'k', ls = ' ')

	#
	Name1, L_V, Re_2d, eRp, eRm, sigma_los, e_sigma, M_half, eM_high, eM_low = sf.get_jdata(jlines)
	Name2, Mst = sf.get_jdata_NEW(jlinesM)
	Mstar = Mst * 1e6
	#plt.plot(L_V, Re_2d * to3D, marker = 'o', markerfacecolor = 'w', markeredgecolor = 'k', ls = ' ', markeredgewidth = 1.5)
	plt.plot(Mstar, Re_2d * to3D, marker = 'o', markerfacecolor = 'w', markeredgecolor = 'k', \
		ls = ' ', markeredgewidth = 1.5, markersize = 10)

	# Redo lines based on M*/pc^2
	MtoL = 1

	Mst1 = 10
	Mst2 = 1e6
	Ib30 = afct.surfB(30)
	Ib35 = afct.surfB(32.5)
	Re1_SDSS = sf.get_Re(Mst1, Ib30, MtoL) * to3D
	Re2_SDSS = sf.get_Re(Mst2, Ib30, MtoL) * to3D
	Re1_LSST = sf.get_Re(Mst1, Ib35, MtoL) * to3D
	Re2_LSST = sf.get_Re(Mst2, Ib35, MtoL) * to3D
	sdssX = lsstX = np.array([Mst1, Mst2])
	sdssY = np.array([Re1_SDSS, Re2_SDSS])
	lsstY = np.array([Re1_LSST, Re2_LSST])

	plt.xlabel(r'$\rm M_{\star} \, (M_{\odot})$', fontsize = fntsz)
	plt.ylabel(r'$\rm R_{1/2} \, (pc)$', fontsize = fntsz)
	plt.loglog()
	plt.plot(sdssX, sdssY, lw = 2, ls = '-', color = 'k')
	plt.plot(lsstX, lsstY, lw = 2, ls = 'dashed', color = 'k')
	plt.xlim(100,2.6*10**7)
	plt.ylim(10,1100*to3D)

	xx = np.array([2, 3, 4, 5, 6, 7])
	xtks = 10 ** xx
	ax.set_xticks(xtks)
	xstr = ['$10^{'+str(x)+'}$' for x in xx]
	ax.set_xticklabels(xstr, fontsize = fntsz)
	yy = np.array([1, 2, 3])
	ytks = 10 ** yy
	ax.set_yticks(ytks)
	ystr = ['$10^{'+str(y)+'}$' for y in yy]
	ax.set_yticklabels(ystr, fontsize = fntsz)
 
	xfake = 1e20
	yfake = 1e20
	xfake1 = 1e21
	yfake1 = 1e21
 
 
 
	plt.plot(xfake, yfake, ls = ' ', marker = 'o', markersize = 10, markeredgecolor = 'k', \
		markerfacecolor = 'w', markeredgewidth = 1.5, label = r'$\rm MW \, Dwarfs$')
	plt.errorbar(xfake1, yfake1, xerr = 100, yerr = 200, ls = ' ', marker = 'h', markersize = 12, color = 'k', \
		label = r'$\rm DES \, Candidates$')
	plt.errorbar(xfake1, yfake1, xerr = 100, yerr = 200, ls = ' ', marker = '.', markersize = 12, color = 'k', \
		label = r'$\rm DES \, Year 2 \, Candidates$')
	plt.plot(xfake, yfake, marker = 's', color = 'k', ls = ' ', markersize = 10, label = r'$\rm Hydra \, II$')

	plt.legend(loc = 4, numpoints = 1, frameon = False, fontsize = fntsz * 0.75)

	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(fntsz)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(fntsz)
	for line in ax.get_xticklines() + ax.get_yticklines():
		line.set_markersize(18)
		line.set_markeredgewidth(3)
	for tick in ax.xaxis.get_minor_ticks():
		tick.label1.set_fontsize(fntsz/2)
	for tick in ax.yaxis.get_minor_ticks():
		tick.label1.set_fontsize(fntsz/2)

	ax.tick_params(which='major',width=axwidth,length=axlength+5, pad = 10)
	ax.tick_params(which='minor',width=axwidth,length=axlength, pad = 10)

 
def pre_plot(Fs = 14, axwidth = 3):
	"""
	
	"""

	plt.rc('axes',linewidth=axwidth)
	fig = plt.figure(1, figsize=(Fs, Fs))
	ax = fig.add_subplot(1,1,1)
	ax = plt.gca()
	return ax

 
def post_plot(ax,xx,yy,fntsz = 32, axwidth = 3, axlength = 12, logx = False, logy = False):
     """
 	
     """
     
     plt.xticks(weight = 'bold')
     plt.yticks(weight = 'bold')
     
     if logx == True:
         xtks = 10. ** xx
         ax.set_xticks(xtks)
         xstr = ['$10^{'+str(x)+'}$' for x in xx]
     else:
         xtks = xx
         ax.set_xticks(xtks)
         xstr = ['$'+str(x) +'$' for x in xx]
     	
     ax.set_xticklabels(xstr, fontsize = fntsz)
     
     if logy == True:
         ytks = 10. ** yy
         ax.set_yticks(ytks)
         ystr = ['$10^{'+str(y)+'}$' for y in yy]
     else:	
         ytks = yy
         ax.set_yticks(ytks)
         ystr = ['$'+str(y)+'$' for y in yy]
     	
     ax.set_yticklabels(ystr, fontsize = fntsz)
	
     for tick in ax.xaxis.get_major_ticks():
         tick.label1.set_fontsize(fntsz)
     for tick in ax.yaxis.get_major_ticks():
         tick.label1.set_fontsize(fntsz)
     for line in ax.get_xticklines() + ax.get_yticklines():
         line.set_markersize(18)
         line.set_markeredgewidth(3)
     for tick in ax.xaxis.get_minor_ticks():
         tick.label1.set_fontsize(fntsz/2)
     for tick in ax.yaxis.get_minor_ticks():
         tick.label1.set_fontsize(fntsz/2)

     ax.tick_params(which='major',width=axwidth,length=axlength+5,pad=15)
     ax.tick_params(which='minor',width=axwidth,length=axlength,pad=15)
 
 
 
def get_edata(elines):
    """
    This program accesses Evan's Kirby 2013 data for the dSphs - their Figure 1
    
    """
    N = np.size(elines)
    Name = []
    LV = np.ndarray(N-1)
    Lerr = np.ndarray(N-1)
    Mstar = np.ndarray(N-1)
    Merr = np.ndarray(N-1)
    FeH = np.ndarray(N-1)
    Ferr = np.ndarray(N-1)
    sig = np.ndarray(N-1)
    serr = np.ndarray(N-1)
    for t in range(N-1):
        Name = np.append(Name, elines[t+1].strip().split()[0])
        LV[t] = np.float(elines[t+1].strip().split()[2])
        Lerr[t] = float(elines[t+1].strip().split()[3])
        Mstar[t] = float(elines[t+1].strip().split()[4])
        Merr[t] = float(elines[t+1].strip().split()[5])
        FeH[t] = float(elines[t+1].strip().split()[6])
        Ferr[t] = float(elines[t+1].strip().split()[7])
        sig[t] = float(elines[t+1].strip().split()[8])
        serr[t] = float(elines[t+1].strip().split()[9])
         
    return Name, LV, Lerr, Mstar, Merr, FeH, Ferr, sig, serr
 
 
def get_M31data(elines):
    """
    This program accesses Vargas 2014 data for the dSphs - their Figure 1
    
    """
    N = np.size(elines)
    Name = []
    Mass = np.ndarray(N-1)
    FeH = np.ndarray(N-1)
    Ferr = np.ndarray(N-1)
    NN = np.ndarray(N-1)
    alphaFe5 = np.ndarray(N-1)
    Ferr05 = np.ndarray(N-1)
    N05 = np.ndarray(N-1)
    alphaFe25 = np.ndarray(N-1)
    Ferr25 = np.ndarray(N-1)
    N25 = np.ndarray(N-1)
    alphaFe2 = np.ndarray(N-1)
    Ferr2 = np.ndarray(N-1)
    N2 = np.ndarray(N-1)
    alphaFe15 = np.ndarray(N-1)
    Ferr15 = np.ndarray(N-1)
    N15 = np.ndarray(N-1)
    
    
    for t in range(N-1):
        Name = np.append(Name, elines[t+1].strip().split()[0])
        Mass[t] = np.float(elines[t+1].strip().split()[1])
        FeH[t] = np.float(elines[t+1].strip().split()[2])
        Ferr[t] = float(elines[t+1].strip().split()[3])
        NN[t] = float(elines[t+1].strip().split()[4])
        alphaFe5[t] = float(elines[t+1].strip().split()[5])
        Ferr05[t] = float(elines[t+1].strip().split()[6])
        N05[t] = float(elines[t+1].strip().split()[7])
        alphaFe25[t] = float(elines[t+1].strip().split()[8])
        Ferr25[t] = float(elines[t+1].strip().split()[9])
        N25[t] = float(elines[t+1].strip().split()[10])
        alphaFe2[t] = float(elines[t+1].strip().split()[11])
        Ferr2[t] = float(elines[t+1].strip().split()[12])
        N2[t] = float(elines[t+1].strip().split()[13])
        alphaFe15[t] = float(elines[t+1].strip().split()[14])
        Ferr15[t] = float(elines[t+1].strip().split()[15])
        N15[t] = float(elines[t+1].strip().split()[16])
         
    return Name, Mass, FeH, Ferr, alphaFe5, Ferr05, alphaFe25, Ferr25, alphaFe2, Ferr2, alphaFe15, Ferr15

def get_Cdata(elines):
    """
    This program accesses Collins 2013 data, asterisk on name means rhalf from Kirby 2014
    No asterisk means McConnachie 2012


    
    """
    N = np.size(elines)
    Name = []
    Mhalf = np.ndarray(N-1)
    Mherr = np.ndarray(N-1)
    Mlerr = np.ndarray(N-1)
    MLhalf = np.ndarray(N-1)
    MLherr = np.ndarray(N-1)
    MLlerr = np.ndarray(N-1)
    Vchalf = np.ndarray(N-1)
    Vcherr = np.ndarray(N-1)
    Vclerr = np.ndarray(N-1)
    rhalf = np.ndarray(N-1)
    rhherr = np.ndarray(N-1)
    rhlerr = np.ndarray(N-1)
    
    
    for t in range(N-1):
        Name = np.append(Name, elines[t+1].strip().split()[0])
        Mhalf[t] = np.float(elines[t+1].strip().split()[1])
        Mherr[t] = np.float(elines[t+1].strip().split()[2])
        Mlerr[t] = float(elines[t+1].strip().split()[3])
        MLhalf[t] = float(elines[t+1].strip().split()[4])
        MLherr[t] = float(elines[t+1].strip().split()[5])
        MLlerr[t] = float(elines[t+1].strip().split()[6])
        Vchalf[t] = float(elines[t+1].strip().split()[7])
        Vcherr[t] = float(elines[t+1].strip().split()[8])
        Vclerr[t] = float(elines[t+1].strip().split()[9])
        rhalf[t] = float(elines[t+1].strip().split()[10])
        rhherr[t] = float(elines[t+1].strip().split()[11])
        rhlerr[t] = float(elines[t+1].strip().split()[12])
         
    return Name, Mhalf, Mherr, Mlerr, MLhalf, MLherr, MLlerr, Vchalf, Vcherr, Vclerr, rhalf, rhherr, rhlerr

def get_MgData(mglines):
    """
    This program accesses Vargas 2014 data for the dSphs - their Figure 1
    
    """
    N = np.size(mglines)
    Name = []
    MgH = np.ndarray(N-1)
    errMg = np.ndarray(N-1)
    FeH = np.ndarray(N-1)
    errFe = np.ndarray(N-1)
 
    
    for t in range(N-1):
        Name = np.append(Name, mglines[t+1].strip().split()[0])
        MgH[t] = np.float(mglines[t+1].strip().split()[1])
        errMg[t] = np.float(mglines[t+1].strip().split()[2])
        FeH[t] = float(mglines[t+1].strip().split()[3])
        errFe[t] = float(mglines[t+1].strip().split()[4])
          
    return Name, MgH, errMg, FeH, errFe



def plot_sfh(stars, runname, ldict, cdict, r_restr = 5, mpart = 250, lw = 2, label = None, norm = True, ls = '-'):
	"""
	Plots cumulative SFH and gas temp hist of stars
	
	"""
	stars.restrict_parts(r_restr)
	arr = stars.arrays.age[stars.booly]
	mpart = stars.arrays.mass

	if norm == True:
		ylog = False
	else:
		ylog = True
		
	cf.plot_cumu(arr, mpart, labels = ldict[label], clr = cdict[label], gt = False, ylog = ylog, norm = norm, lw = lw, ls = ls)
	plt.legend(loc = 4, numpoints = 1, frameon = False, fontsize = 16)


def plot_gasT_hist(gas, runname, ldict, cdict, r_restr = 5, lw = 2, cumu = False):
	"""
	Plots gas temp hist of stars
	
	"""
# 	gas.restrict_parts(r_restr)
	
	if cumu == False:
		plt.hist(np.log10(gas.arrays.temp), histtype = 'step', normed = True, lw = 2, bins = 50, color = cdict['%s'%(runname)], label = ldict['%s'%(runname)])
		plt.xlabel(r"$\rm log \, T \, (K)$", fontsize = 16)
		plt.ylabel(r"$\rm N$", fontsize = 16)
		plt.legend(loc = 1, numpoints = 1, frameon = False, fontsize = 12)
		plt.xlim(2,5)
	else:
		cf.plot_cumu(np.log10(gas.arrays.temp), mass = 250, labels = ldict['%s'%(runname)], clr = cdict['%s'%(runname)], loc = 4, gt = False, norm = True)
		plt.xlabel(r"$\rm log \, T \, (K)$", fontsize = 16)
		plt.ylabel(r"$\rm f(>T)$", fontsize = 16)
		plt.legend(loc = 4, numpoints = 1, frameon = False, fontsize = 12)
		plt.xlim(2,5)
	
################## Metallicity Stuff #######################
def plot_alpha(st_gs, r_restr, alph = 'Ca', vlims = [100,1e5], xlims = [-3.5,1], ylims = [-0.6,0.4]):
	"""
	Plots alpha (currently Ca, change  with 'alph') of gas or stars
	
	"""
	st_gs.restrict_parts(r_restr)
	st_gs.plot_alpha_Fe_vs_Fe_H_hist(lims = vlims, alph = alph)
	plt.xlim(xlims)
	plt.ylim(ylims)



################## Density-T #######################
def plot_dens_T(gas, r_restr = 5., vlims = [100,1e5], xlims = None, ylims = None):
# 	gas.restrict_parts(r_restr)
	gas.plot_T_vs_density_hist(lims = vlims)
	
	if xlims is not None:
		plt.xlim(xlims)
		plt.ylim(ylims)



################## SFR / Variance ##################################
def plot_sfr(stars, runname, cdict, r_restr=5., smooth1=10, smooth2 = 200, clr = 'k'):
	"""
	This function plots the SFR
	
	"""
	
	stars.restrict_parts(5.)
	cf.plot_sfr(stars, smooth = smooth2, clr = cdict['%s'%(runname1)], specific = False, lw = 2)

def plot_sfr_fburst(stars, runname, cdict, r_restr=5., smooth1=10, smooth2 = 200):
	"""
	Plots the fraction of sfr in bursts, postbursts and in between
	
	"""
	fburst = cf.compare_sfr(stars, smooth1 = smooth1, smooth2 = smooth2, fac = fac, clr = cdict['%s'%(runname)])
	plt.annotate("f_burst = %.2f"%(fburst), xy = (13.5,-1.2), fontsize = 12, color = cdict['%s'%(runname)])


################## Density Profiles with Power radius ####################
def plot_dens_power(dm, runname, label = None, r_restr = 500, nbins = 100, alpha = 0.2):

	"""
	Plots density with power radius (200 DM particles)

	"""
	dm.orderp()
	power = dm.arrays.r[199]
	dm.reset_to_orig()
	dm.recenter()
	dm.restrict_parts(r_restr)
	dm.plot_rho(nbins = nbins, logbins = True, logR = True, clr = cdict['%s'%(runname)], lw = 2.5, labels = label)
	plt.axvspan(-3, np.log10(power), facecolor=cdict['%s'%(runname)], edgecolor = cdict['%s'%(runname)], alpha=alpha)
	

################## 2D Histograms #######################

def plot_2d_hist(dm_or_s, r_restr = 60, r_ext = 60, h_ext = 60, dmbins = 1500, sbins = 800, \
	xlims = [-30,30], ylims = [-30,30], cmap = 'bone', topdown = True, wt = 'mass', lims = [-1,-1]):
	"""
	2D hist of either dm or stars
	
	"""
	if r_restr is not None:
		dm_or_s.restrict_parts(r_restr)
	
	if wt == 'mass':
		dm_or_s.plot_mass_2Dhist(r_ext = r_restr, h_ext = h_ext, bins = dmbins, cmap = cmap, topdown = topdown)
	else:
		dm_or_s.plot_alpha_2Dhist(r_ext = r_restr, h_ext = h_ext, bins = sbins, cmap = cmap, topdown = topdown, lims = lims)
	plt.xlim(xlims)
	plt.ylim(xlims)
	

def plot_2d_hist_basic(dm, r_restr = None, bins = 500, topdown = True, cmap = 'bone', cmin = 0, xylims = None):
	"""
	2D hist of dm only simple
	
	"""
	if r_restr is not None:
		dm.restrict_parts(r_restr)
	fig = plt.figure(1, figsize = (10,10))
	ax1 = fig.add_subplot(1,1,1)
	ax1.patch.set_facecolor('black')
	x = dm.arrays.x[dm.booly]
	y = dm.arrays.y[dm.booly]
	z = dm.arrays.z[dm.booly]
	weights = dm.arrays.mass[dm.booly]
	if topdown == True:
		ax1.hist2d(x, y, norm=matplotlib.colors.LogNorm(), bins = bins, weights = weights, cmap = cmap, cmin = cmin)
	else:
		ax1.hist2d(x, z, norm=matplotlib.colors.LogNorm(), bins = bins, weights = weights, cmap = cmap, cmin = cmin)
	if xylims is not None:
		plt.xlim(xylims[0],xylims[1])
		plt.ylim(xylims[2],xylims[3])
	plt.show()

	return ax1

#########  MDF #########
def plot_mdf(stars, runname, cdict, label, r_restr = 5, bins = 50):
	"""
	Metallicity distribution function
	
	"""
	
	stars.restrict_parts(r_restr)
	plt.hist(stars.Fe_H[stars.booly], histtype = 'step', normed = True, lw = 2, bins = bins, color = cdict['%s'%(runname)], label = label)
	plt.legend(loc = 1, numpoints = 1, frameon = False, fontsize = 16)

#########  MDF #########
def plot_Fe_grad(stars, runname, cdict, label, r_restr = 5, bins = 50):
	"""
	Metallicity distribution function
	
	"""
	stars.Fe_grad(D = r_restr, Fs = 8, clr = cdict['%s'%(runname)])
	plt.legend(loc = 1, numpoints = 1, frameon = False, fontsize = 16)

def open(fname):
	"""
	Uses pickle to load object
	
	"""
	with open("sA_m0930_pickle", 'rb') as output: 
            return pickle.load(output)

def get_flat_angle(halo, plane='xy', nbins = 30, center = True):
    """
    This function rotates glx until ellipticity minimized then that is the angle for rotation to calc
    Ellipticity
    """

    mass = halo.arrays.mass[halo.booly]
    pos = halo.pos[halo.booly]
    vel = halo.vel[halo.booly]
    
    X = halo.rvir/10.
        
    # Now calculate angular momentum vector of galaxy using sum of cross products of all points
    ##############
    Tratio = np.ndarray(nbins)
    q = np.ndarray(nbins)
    N = len(mass)
       
    alltheta = np.linspace(0,np.pi,nbins)
    for i in range(nbins):
	    rot_theta = alltheta[i]
        # Now get los velocity for given rotation angle and plane
	    if plane == 'xy':
		    vlos = vel[:,2]
		    wd = pos[:,0] * np.cos(rot_theta) - pos[:,1] * np.sin(rot_theta)
		    ht = pos[:,1] * np.cos(rot_theta) + pos[:,0] * np.sin(rot_theta)
	    elif plane == 'xz':
		    vlos = -vel[:,1]
		    wd = pos[:,0] * np.cos(rot_theta) - pos[:,2] * np.sin(rot_theta)
		    ht = pos[:,2] *  np.cos(rot_theta) + pos[:,0] * np.sin(rot_theta)
	    elif plane == 'yz':
		    vlos = vel[:,0]
		    wd = pos[:,1] * np.cos(rot_theta) - pos[:,2] * np.sin(rot_theta)
		    ht = pos[:,2] * np.cos(rot_theta) + pos[:,1] * np.sin(rot_theta)
    
    
	    wbins = np.linspace(0,X,nbins+1)
	    hbins = np.linspace(0,X,nbins+1)
	    
	    vrot = np.ndarray((nbins,nbins))
	    sig = np.ndarray((nbins,nbins))
	    ninbin = np.ndarray((nbins,nbins))
	    ratio = np.ndarray((nbins,nbins))
	    xmidbin = np.ndarray((nbins,nbins))
	    ymidbin = np.ndarray((nbins,nbins))

         
	    for t in range(nbins):
		    for y in range(nbins):
			    ninbin[t][y] = len(vlos[(((wd >= wbins[y])&(wd < wbins[y+1]))&((ht >= hbins[t])&(ht < hbins[t+1])))])
			    xmidbin[t][y] = (wbins[y] + wbins[y+1]) / 2.0
			    ymidbin[t][y] = (hbins[t] + hbins[t+1]) / 2.0
			    if ninbin[t][y] > 5:
				    vrot[t][y] = np.mean(vlos[(((wd >= wbins[y])&(wd < wbins[y+1]))&((ht >= hbins[t])&(ht < hbins[t+1])))])
				    sig[t][y] = np.std(vlos[(((wd >= wbins[y])&(wd < wbins[y+1]))&((ht >= hbins[t])&(ht < hbins[t+1])))])
				    if sig[t][y] > 0:
					    ratio[t][y] = vrot[t][y] / sig[t][y]
				    else:
					    ratio[t][y] = vrot[t][y]
			    else:
				    vrot[t][y] = sig[t][y] = ratio[t][y] = 0

	    Tratio[i], q[i] = sf.get_sauron_ratio_q(ninbin, vrot, sig, ratio, xmidbin, ymidbin)   

    return Tratio[np.argmin(q)], 1 - np.amin(q), alltheta[np.argmin(q)], N


def get_obsdat(lines):
    """
    obsdat has no delim
    
    """
    N = len(lines[0].strip().split())
    M = 40 # M-1 dwarfs
    cols = np.ndarray((M,N))
    colnames = lines[0].strip().split()
    dwarfnames = np.ndarray(M, dtype = '|S32')
    dwarfnames[0] = lines[0].strip().split()[0]
    for m in range(M-1):
        m += 1
        dwarfnames[m] = lines[m].strip().split()[0]
        stringarray = lines[m].strip().split()[1:]
        cols[m][1:] = np.array([float(y) for y in stringarray])
    
    obsdat = rec.array([dwarfnames[1:], cols[1:,1], cols[1:,2], cols[1:,3], cols[1:,4], cols[1:,5],\
        cols[1:,6], cols[1:,7], cols[1:,8], cols[1:,9], cols[1:,10],cols[1:,11], cols[1:,12], cols[1:,13], \
        cols[1:,14], cols[1:,15], cols[1:,16], cols[1:,17], cols[1:,18], cols[1:,19], cols[1:,20], \
        cols[1:,21], cols[1:,22], cols[1:,23], cols[1:,24], cols[1:,25], cols[1:,26]], \
        dtype = [('Name', '|S32'), ('id', '<i8'), \
        ('rmed', '<f8'),('rmerr_m', '<f8'), ('rmerr_p', '<f8'), ('sig', '<f8'),  \
        ('serr_m', '<f8'), ('serr_p', '<f8'), ('vrot', '<f8'), ('verr_m', '<f8'), \
        ('verr_p', '<f8'), ('posang', '<f8'), ('aerr_m', '<f8'), ('aerr_p', '<f8'), \
        ('ellip', '<f8'), ('eerr_m', '<f8'), ('eerr_p', '<f8'), ('mass', '<f8'), \
        ('dist', '<f8'), ('nstar', '<i8'), ('bayes', '<f8'), ('ratio', '<f8'), ('rerr_m', '<f8'),\
        ('rerr_p', '<f8'), ('rhalf_full', '<f8'), ('rhalf_m', '<f8'), ('rhalf_p', '<f8')])
    
    return obsdat



def get_rpisodat(fname, delim = ',', Num = 40):
    """
    Shortcut to get any file with specified delimiter
    
    """
    
    dat = open(fname, 'r')
    lines = dat.readlines()
    N = len(lines[0].strip(',\n').split(delim))
    M = Num # M-1 dwarfs
    cols = np.ndarray((M,N))
    colnames = lines[0].strip(',\n').split(delim)
    dwarfnames = np.ndarray(M, dtype = '|S32')
    dwarfnames[0] = lines[0].strip(',\n').split(delim)[0]
    for m in range(M-1):
        m += 1
        dwarfnames[m] = lines[m].strip(',\n').split(delim)[0]
        stringarray = lines[m].strip(',\n').split(delim)[1:]
        cols[m][1:] = np.array([float(y) for y in stringarray])
    
        
    rpisodat = rec.array([dwarfnames[1:], cols[1:,1], cols[1:,2], cols[1:,3], cols[1:,4], cols[1:,5],\
        cols[1:,6], cols[1:,7], cols[1:,8], cols[1:,9], cols[1:,10],cols[1:,11], cols[1:,12], cols[1:,13], \
        cols[1:,14], cols[1:,15], cols[1:,16], cols[1:,17], cols[1:,18], cols[1:,19], cols[1:,20], \
        cols[1:,21], cols[1:,22], cols[1:,23], cols[1:,24], cols[1:,25], cols[1:,26], cols[1:,27], cols[1:,28], \
        cols[1:,29], cols[1:,30], cols[1:,31], cols[1:,32], cols[1:,33], cols[1:,34], cols[1:,35], cols[1:,36], cols[1:,37], \
        cols[1:,38], cols[1:,39], cols[1:,40], cols[1:,41], cols[1:,42], cols[1:,43], cols[1:,44]], \
        dtype = [('Name', '|S32'), ('rhalf', '<f8'),('ellip', '<f8'), ('rhalf_sph', '<f8'), \
        ('r75', '<f8'), ('r90', '<f8'), ('Bayes', '<f8'), ('ratio_flat', '<f8'), ('flat_m', '<f8'), ('flat_p', '<f8'), \
        ('ratio_rhalf', '<f8'), ('rhalf_m', '<f8'), ('rhalf_p', '<f8'), ('ratio_rhalf_sph', '<f8'), \
        ('sph_m', '<f8'), ('sph_p', '<f8'), ('ratio_r75', '<f8'), ('r75_m', '<f8'), \
        ('r75_p', '<f8'), ('ratio_r90', '<f8'), ('r90_m', '<f8'), ('r90_p', '<f8'), \
        ('sigma_vrot', '<f8'), ('svr_m', '<f8'), ('svr_p', '<f8'), ('vrot_rhalf', '<f8'), ('vrrh_m', '<f8'), \
        ('vrrh_p', '<f8'), ('vrot_sph', '<f8'), ('vrsph_m', '<f8'), ('vrsph_p', '<f8'), ('vrot_r75', '<f8'), \
        ('vrr75_m', '<f8'), ('vrr75_p', '<f8'), ('vrot_r90', '<f8'), ('vrr90_m', '<f8'), ('vrr90_p', '<f8'), \
        ('vrot_flat', '<f8'), ('vrf_m', '<f8'), ('vrf_p', '<f8'), ('sigma_flat', '<f8'), ('sf_m', '<f8'), \
        ('sf_p', '<f8'), ('Bayes_flat', '<f8'), ('Bayes_PISO', '<f8')])
    
    return rpisodat
    

def get_rvardata(elines):
    """
   R_75  R_90  R_95  v_sigma_75  v_sigma_75_minus  v_sigma_75_plus  v_sigma_90  
   v_sigma_90_minus  v_sigma_90_plus  v_sigma_95  v_sigma_95_minus  v_sigma_95_plus  Bayes  
    """
    N = np.size(elines)
    Name = np.ndarray(N-1, dtype = '|S32')
    R_75 = np.ndarray(N-1)
    R_90 = np.ndarray(N-1)
    R_95 = np.ndarray(N-1)
    v_sigma_75 = np.ndarray(N-1)
    v_sigma_75_minus = np.ndarray(N-1)
    v_sigma_75_plus = np.ndarray(N-1)
    v_sigma_90 = np.ndarray(N-1)
    v_sigma_90_minus = np.ndarray(N-1)
    v_sigma_90_plus = np.ndarray(N-1)
    v_sigma_95 = np.ndarray(N-1)
    v_sigma_95_minus = np.ndarray(N-1)
    v_sigma_95_plus = np.ndarray(N-1)
    Bayes_rvar = np.ndarray(N-1)
    for t in range(N-1):
        Name[t] = elines[t+1].strip().split()[0]
        R_75[t] = np.float(elines[t+1].strip().split()[1])
        R_90[t] = np.float(elines[t+1].strip().split()[2])
        R_95[t] = np.float(elines[t+1].strip().split()[3])
        v_sigma_75[t] = np.float(elines[t+1].strip().split()[4])
        v_sigma_75_minus[t] = np.float(elines[t+1].strip().split()[5])
        v_sigma_75_plus[t] = np.float(elines[t+1].strip().split()[6])
        v_sigma_90[t] = np.float(elines[t+1].strip().split()[7])
        v_sigma_90_minus[t] = np.float(elines[t+1].strip().split()[8])
        v_sigma_90_plus[t] = np.float(elines[t+1].strip().split()[9])
        v_sigma_95[t] = np.float(elines[t+1].strip().split()[10])
        v_sigma_95_minus[t] = np.float(elines[t+1].strip().split()[11])
        v_sigma_95_plus[t] = np.float(elines[t+1].strip().split()[12])
        Bayes_rvar[t] = np.float(elines[t+1].strip().split()[13])

        
        alldat = rec.array([Name, R_75, R_90, R_95, v_sigma_75, v_sigma_75_minus, \
            v_sigma_75_plus, v_sigma_90, v_sigma_90_minus, v_sigma_90_plus, v_sigma_95, \
            v_sigma_95_minus, v_sigma_95_plus, Bayes_rvar], dtype = [('Name', '|S32'), \
            ('R75', '<f8'),('R90', '<f8'), ('R95', '<f8'), ('v_sigma_75', '<f8'), ('v_sigma_75_minus', '<f8'), \
            ('v_sigma_75_plus', '<f8'), ('v_sigma_90', '<f8'), ('v_sigma_90_minus', '<f8'), \
            ('v_sigma_90_plus', '<f8'), ('v_sigma_95', '<f8'), ('v_sigma_95_minus', '<f8'), \
            ('v_sigma_95_plus', '<f8'), ('Bayes_rvar', '<f8')])

    return alldat




def plot_sfr(s, ax1, smooth, clr = None, specific = True, lw = 1, ls = '-', axes = True, \
	Fs = 12, axwidth = 3, axlength = 4, fntsz = 10):
	"""
	Plots the SFR smoothed over timescale smooth
	
	"""
	
	ax1 = plt.gca()
	plt.rc('axes',linewidth=axwidth)
	
	from astropy.cosmology import FlatLambdaCDM
	cosmo = FlatLambdaCDM(Om0=s.header['omega_matter'], H0=s.header['hubble']*100)
	
	if clr is None:
		clr = s.color

	ttot = 13.8
	nbins = np.int(np.floor(ttot / (smooth/1000.)))
	bins = np.linspace(0,13.8,nbins)
	binavg_or_sum, midbins, varbins, rbins = cf.bin_data(s.arrays.age[s.booly], \
														  s.arrays.mass[s.booly], nbins, logbins = False, suminbin = True)
	peryr = binavg_or_sum / (ttot*1e9/nbins)
	mtot = np.ndarray(nbins)
	mtot[0] = 0# np.sum(disk.arrays.mass)
	for i in range(nbins - 1):
		t = i+1
		mtot[t] = mtot[t-1] + binavg_or_sum[t]
	if smooth < 1000:
		label = r"$\rm %i \, Myr \, Avg$"%(smooth)
	else:
		label = r"$\rm %i \, Gyr \, Avg$"%(smooth/1000)
	if specific == True:
		logperyr = np.log10(peryr/mtot)
	else:
		logperyr = np.log10(peryr)
	logperyr[logperyr == -np.inf] = -10
	if specific == True:
		plt.plot(midbins, logperyr[::-1], label = label, color = clr, lw = lw, ls = ls)
		plt.ylabel(r"$\rm log \,sSFR \, [yr^{-1}]$", fontsize = fntsz+6)
	else:
		plt.plot(midbins, logperyr[::-1], label = label, color = clr, lw = lw, ls = ls)
		plt.ylabel(r"$\rm log \, SFR \, [M_* \, / \, yr^{-1}]$", fontsize = fntsz+6)
		
	plt.legend(loc = 1, numpoints = 1, frameon = False, fontsize = fntsz+6)
	plt.ylim(-5, -2)
	plt.xlim(13.8,0)
	
	if axes == True:
	
		zs = np.array([10,5,3,2,1,0.5,0.25,0])
		aa = [cosmo.lookback_time(z).value for z in zs]
		
		ax2 = ax1.twiny()
		ax2.set_xticks(aa)
		ax2.set_xticklabels(['{:g}'.format(z) for z in zs])
		ax2.set_xlim(ax1.get_xlim())
		zstr = ['$'+str(z)+'$' for z in zs]
		ax2.set_xticklabels(zstr, fontsize=fntsz)
		ax2.set_xlabel(r'$\rm z$', fontsize = fntsz)
	
		xx = [12,10,8,6,4,2,0]
		xstr = ['$'+str(x)+'$' for x in xx]
		ax1.set_xticklabels(xstr[::-1], fontsize = fntsz)
		
		yy = [-5.0,-4.5,-4.0,-3.5,-3.0,-2.5,-2.0]
		ystr = ['$'+str(y)+'$' for y in yy]
		
		ax1.set_yticklabels(ystr, fontsize = fntsz)
		plt.xticks(weight = 'bold')
		plt.yticks(weight = 'bold')
		
		
		ax1.set_ylabel(r'$\rm log \, SFR \, (M_{\odot} yr^{-1})$', fontsize = fntsz)
		ax1.set_xlabel(r"$\rm Lookback \, Time \, (Gyr)$", fontsize = fntsz)
		
		ax = ax1
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fntsz)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fntsz-2)
		for line in ax.get_xticklines() + ax.get_yticklines():
			line.set_markersize(fntsz)
			line.set_markeredgewidth(axwidth)
		for tick in ax.xaxis.get_minor_ticks():
			tick.label1.set_fontsize(fntsz/2)
		for tick in ax.yaxis.get_minor_ticks():
			tick.label1.set_fontsize(fntsz/2)
			
		ax.tick_params(which='major',width=axwidth,length=axlength+2,pad = 10)
		ax.tick_params(which='minor',width=axwidth,length=axlength,pad=10)

		#####################

		ax = ax2
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fntsz)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fntsz-2)
		for line in ax.get_xticklines() + ax.get_yticklines():
			line.set_markersize(fntsz)
			line.set_markeredgewidth(axwidth)
		for tick in ax.xaxis.get_minor_ticks():
			tick.label1.set_fontsize(fntsz/2)
		for tick in ax.yaxis.get_minor_ticks():
			tick.label1.set_fontsize(fntsz/2)
			
		ax.tick_params(which='major',width=axwidth,length=axlength+2,pad=10)
		ax.tick_params(which='minor',width=axwidth,length=axlength,pad=10)


def plot_rho(r,m, nbins=100, logbins = False, logR=False, clr = 'k', lw = 2, ls = '-', labels = None):
	PI = np.pi
	mpart = m[np.argsort(r)] # Get everything in order of increasing r
	rad = r[np.argsort(r)]
	power = rad[1999]
	rbins = gen_bins(rad, nbins, logbins)

	vperbin = np.ndarray(nbins) # will hold volume per bin
	mperbin = np.ndarray(nbins) # will hold mass per bin
	nperbin = np.ndarray(nbins) # will hold number of particles per bin
	midbin = np.ndarray(nbins) # will store middle of r bin

	# Calculate volume per bin
	for c in range(nbins):
		vperbin[c] = (rbins[c+1] ** 3 - rbins[c] ** 3) * (4./3) * PI # volume per bin (kpc^3)
		radinbin = rad[((rad >= rbins[c]) & (rad < rbins[c+1]))]
		massinbin = mpart[((rad >= rbins[c]) & (rad < rbins[c+1]))] # nec for diff masses
		nperbin[c] = np.size(radinbin) # number of particles in each bin
		mperbin[c] = np.sum(massinbin)
		midbin[c] = (rbins[c+1] + rbins[c]) / 2

	# Volume density = npart * mpart / vbin
	rho = mperbin / vperbin    
	if logR == True:
		plt.plot(midbin, rho, color = clr, linewidth = lw, linestyle = ls, label = labels) # Pass string for legend
		plt.xlabel(r'$\rm R (kpc)$', fontsize = 22)
		plt.loglog()
	else:
		plt.plot(midbin, np.log10(rho), color = clr, linewidth = lw, linestyle = ls, label = labels) # Pass string for legend
		plt.xlabel(r'$\rm R (kpc)$', fontsize = 22)
	return midbin,rho,mpart,rad, power

