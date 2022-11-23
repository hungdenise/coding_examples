# all sextractor fits
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
import numpy as np
import math
from scipy import integrate
from scipy.optimize import curve_fit
import os

plt.rcParams.update({'font.size': 20})
plot = True
links = ['1.5', '2', '2.5', '3']
fields = ['COSMOS', 'CFHTLS-D1', 'ECDFS']
field_err = {'ECDFS': 0.034747673, 'CFHTLS-D1': 0.035024802, 'COSMOS': 0.032169747} # average error in field from clipped sigma(delta_z/(1+z))
detects = ['3', '3.5']
mareas = ['100']

def gaus(x,a,x0,sigma):
  "1-d gaussian: gaus(x, a, x0, sigma)"
  return a*np.exp(-(x-x0)**2/(2*sigma**2))

cosmo = lambda z: 1.0/math.sqrt(omega_m*(1+z)**3 + omega_lambda)
H0 = 70.0
omega_m = 0.27
omega_lambda = 0.73

pos = 'known_pos.txt'
pdt = np.dtype({'names':['field', 'struct_id', 'ra', 'dec', 'z1', 'z2', 'delta_prop_d', 'v', 'num_redshifts', 'radius'], 'formats':['U11',np.int,np.float,np.float,np.float,np.float,np.float,np.float,np.int,np.float]})
post = np.loadtxt(pos,dtype=pdt)

pfield = post['field']
pname = post['struct_id']
pra = np.deg2rad(post['ra'])
pdec = np.deg2rad(post['dec'])
pz1 = post['z1']
pz2 = post['z2']

pz = (pz1 + pz2) * 0.5

spec_bandcut = 'bandcut_spec.txt'
phot_bandcut = 'bandcut_phot.txt'
sdt = np.dtype({'names':['field', 'ra', 'dec', 'z', 'pid'], 'formats':['U10', np.float,np.float,np.float,np.int]})
specs_cut = np.loadtxt(spec_bandcut,dtype=sdt)
phot_cut = np.loadtxt(phot_bandcut,dtype=sdt)

spec_field = specs_cut['field']
spec_ras = specs_cut['ra']
spec_decs = specs_cut['dec']
spec_zs = specs_cut['z']
sids = specs_cut['pid']

phot_field = phot_cut['field']
phot_ras = phot_cut['ra']
phot_decs = phot_cut['dec']
phot_zs = phot_cut['z']
pids = phot_cut['pid']

for lk in links:
  if '.' in lk:
    link = lk + '0000'
  else:
    link = lk
  for field in fields:
    fcut = np.where(spec_field == field)
    spec_ra = np.deg2rad(spec_ras[fcut])
    spec_dec = np.deg2rad(spec_decs[fcut])
    spec_z = spec_zs[fcut]
    sid = sids[fcut]
    
    fcut = np.where(phot_field == field)
    phot_ra = np.deg2rad(phot_ras[fcut])
    phot_dec = np.deg2rad(phot_decs[fcut])
    phot_z = phot_zs[fcut]
    pid = pids[fcut]
    
    plt.figure(figsize=(12, 12))
    plt.scatter(phot_ra, phot_dec, s=30, label='z-phot')
    plt.scatter(spec_ra, spec_dec, s=30, label='z-spec')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.legend(loc='best')
    plt.savefig('spec_phot_' + field + '.png', bbox_inches='tight')
    
    for detect in detects:
      for marea in mareas:
        f = open('regmatch/regmatch_' + field + '_d' + detect + '_a' + marea + '_db01_' + link + 'd.txt', 'r')
        lines = f.readlines()
        f.close()
        
        no = []
        detect_no = []
        z = []
        flux = []
        flux_err = []
        ra = []
        dec = []
        area = []
        n = -1
        
        edge_frac = []
        edge = 0
        for l in lines[1:]:
          row = l.split(', ')
          if len(row) > 2:
            no[n].append('f' + row[0] + 's' + row[1])
            detect_no[n].append(row[1])
            z[n].append(row[2])
            flux[n].append(row[4])
            flux_err[n].append(row[5])
            ra[n].append(row[6])
            dec[n].append(row[7])
            area[n].append(row[3])
            if '-' in row[1]:
              edge += 1
          else:
            if n > -1:
              edge_frac.append(edge / len(no[n]))
              edge = 0
            n += 1
            no.append([])
            detect_no.append([])
            z.append([])
            flux.append([])
            flux_err.append([])
            ra.append([])
            dec.append([])
            area.append([])
  
        subsets = set([])
        # exclude cases with detections with very large flux errors or outside search limits
        for i in range(n):
          struct_ra = np.asarray(ra[i], dtype=np.float64)
          struct_dec = np.asarray(dec[i], dtype=np.float64)
          struct_flux_err = np.asarray(flux_err[i], dtype=np.float64)
          large_err = np.where(struct_flux_err > 1000)
          if len(struct_flux_err[large_err]) > 0: 
            subsets.add(i)
        		  
        # look for structure subsets and exclude those from fitting to avoid duplicates
        for i in range(n):
          if len(no[i]) < 5:
            subsets.add(i)
          elif i not in subsets:
            set1 = set([])
            no_sub = no[i]
            z_sub = z[i]
            for x in range(len(no[i])):
              pair1 = (no_sub[x],z_sub[x])
              set1.add(pair1)
            for j in range(i+1,n):
              set2 = set([])
              no_sub2 = no[j]
              z_sub2 = z[j]
              for x in range(len(no[j])):
                pair2 = (no_sub2[x],z_sub2[x])
                set2.add(pair2)
              if set2.issubset(set1):
                subsets.add(j)
        	
        sigma = 0.01
        
        passed_ids = []
        r2s = [] # keep track of goodness of fit
        
        total_passed = 0
        total_links = 0
        for i in range(n):
          # only fit for explicit non-subset structures with at least 5 points
          if i not in subsets and edge_frac[i] < 1:
            struct_z = np.array(z[i], dtype=np.float64)
            struct_area = np.array(area[i], dtype=np.int)
            struct_flux = np.array(flux[i], dtype=np.float64)
            struct_flux_err = np.array(flux_err[i], dtype=np.float64)
            struct_ra = np.array(ra[i], dtype=np.float64)
            struct_dec = np.array(dec[i], dtype=np.float64)
        
            flux_area = struct_flux * struct_area
        
            mean_ra = np.sum(struct_flux * struct_ra) / np.sum(struct_flux)
            mean_dec = np.sum(struct_flux * struct_dec) / np.sum(struct_flux)
            
            a = np.amax(struct_flux)
            z_peak = struct_z[np.where(struct_flux == a)][0] # redshift of peak value of cluster
        
            # do the gaussian fitting
            try:
              popt,pcov = curve_fit(gaus,struct_z,struct_flux,p0=[a,z_peak,sigma],sigma=struct_flux_err,absolute_sigma=True)
              total_links += 1
              
              mean = popt[1]
              std = np.abs(popt[2])
              
              gauss_amp = popt[0]
              gauss_area = gauss_amp * np.sqrt(np.pi * 2.0 * std**2)
        
              perr = np.sqrt(np.diag(pcov))
              gauss_amp_err = perr[0]
              mean_err = perr[1]
              std_err = perr[2]
              errs = gauss_area * np.sqrt((gauss_amp_err/gauss_amp)**2 + (std_err/std)**2)
              if gauss_amp < np.max(struct_flux) * 1.2 and np.abs(mean - np.median(struct_z)) < 0.02 and std < 0.06 and mean_err < 0.02:
                # residual sum of squares
                ss_res = np.sum((struct_flux - gaus(struct_z,gauss_amp,mean,std)) ** 2)
                
                # total sum of squares
                ss_tot = np.sum((struct_flux - np.mean(struct_flux)) ** 2)
                
                # r-squared
                r2 = 1 - (ss_res / ss_tot)
            
                passed_ids.append(i)
                r2s.append(r2)
            except RuntimeError:
              continue
            except TypeError:
              print('Unexpected error in input, check:')
              print(i)
              print(type(struct_z), struct_z)
              print(type(struct_flux), struct_flux)
              print(type(struct_flux_err), struct_flux_err)
              print(a,z_peak,sigma)
              exit()
            except ValueError:
              continue
        
        linked = []
        
        test_r2s = r2s[:]
        test_ids = passed_ids[:]
        while r2s:
          max_r2 = np.max(r2s)
          max = np.where(test_r2s == max_r2)[0][0]
        
          this_i = test_ids[max]
          sextractor_obj = no[this_i]
        
          unique_link = True
          for s in sextractor_obj:
            if s in linked:
              passed_ids.remove(this_i)  # removes the first item from the list whose value is this_i
              unique_link = False
              break # stop search once match is found
        
          r2s.remove(max_r2) # removes the first item from the list whose value is max_r2
          if unique_link:
            linked += sextractor_obj
        
        total_passed += len(passed_ids)
        #print(field, detect, marea, len(passed_ids))
        for i in passed_ids:
          # only fit for explicit non-subset structures with at least 5 points
          this_edge_frac = edge_frac[i]
          struct_no = np.array(detect_no[i], dtype=np.int)
          struct_z = np.asarray(z[i], dtype=np.float64)
          struct_area = np.asarray(area[i], dtype=np.int)
          struct_flux = np.asarray(flux[i], dtype=np.float64)
          struct_flux_err = np.asarray(flux_err[i], dtype=np.float64)
          struct_ra = np.asarray(ra[i], dtype=np.float64)
          struct_dec = np.asarray(dec[i], dtype=np.float64)
        
          flux_area = struct_flux * struct_area
        
          mean_ra = np.sum(struct_flux * struct_ra) / np.sum(struct_flux)
          mean_dec = np.sum(struct_flux * struct_dec) / np.sum(struct_flux)
          
          a = np.amax(struct_flux)
          z_peak = struct_z[np.where(struct_flux == a)][0] # redshift of peak value of cluster
          sigma = 0.01
        
          # do the gaussian fitting
          popt,pcov = curve_fit(gaus,struct_z,struct_flux,p0=[a,z_peak,sigma],sigma=struct_flux_err,absolute_sigma=True)
          
          mean = popt[1]
          std = np.abs(popt[2])
          
          gauss_amp = popt[0]
          gauss_area = gauss_amp * np.sqrt(np.pi * 2.0 * std**2)
        
          perr = np.sqrt(np.diag(pcov))
          gauss_amp_err = perr[0]
          mean_err = perr[1]
          std_err = perr[2]
          errs = gauss_area * np.sqrt((gauss_amp_err/gauss_amp)**2 + (std_err/std)**2)
          
          directory = 'fits/nofit_' + field + '_c80_d' + detect + '_a' + marea + '_db01_link' + lk + '/'
          if not os.path.exists(directory):
            os.makedirs(directory)
          Mpc = 299792./H0*integrate.quad(cosmo, 0, mean)[0]
          dA = Mpc/(1+mean)
          radius_Mpc = np.rad2deg(1.0/dA)
        
          d_spec = np.zeros(len(spec_ra))
          d_phot = np.zeros(len(phot_ra))
          
          for x in range(len(spec_ra)):
            d_spec[x] = math.acos(math.sin(spec_dec[x]) * math.sin(np.deg2rad(mean_dec)) + math.cos(spec_dec[x]) * math.cos(np.deg2rad(mean_dec)) * math.cos(spec_ra[x]-np.deg2rad(mean_ra)))
          for x in range(len(phot_ra)):
            d_phot[x] = math.acos(math.sin(phot_dec[x]) * math.sin(np.deg2rad(mean_dec)) + math.cos(phot_dec[x]) * math.cos(np.deg2rad(mean_dec)) * math.cos(phot_ra[x]-np.deg2rad(mean_ra)))
          
          d_spec = d_spec * dA
          d_phot = d_phot * dA
        
          # 2 Mpc for protoclusters
          gals_spec = np.where(d_spec < 2.0)
          gals_phot = np.where(d_phot < 2.0)
          
          phot_z1 = phot_z[gals_phot]
          spec_z1 = spec_z[gals_spec]
          pid1 = pid[gals_phot]
          sid1 = sid[gals_spec]
          
          delta_redshift = field_err[field] * (1 + mean) # average error * (1+z)
          selects = np.where((spec_z1 > mean - delta_redshift) & (spec_z1 < mean + delta_redshift))
          selectp = np.where((phot_z1 > mean - delta_redshift) & (phot_z1 < mean + delta_redshift))
          selects2 = np.where((spec_z1 <= mean - delta_redshift) | (spec_z1 >= mean + delta_redshift))
          match_inrange = 0
          bad_phot_match = 0
          
          sid2 = sid1[selects]
          pid2 = pid1[selectp]
          sid3 = sid1[selects2]
          pid_nomatch = []
          
          flag = 'nomatch'
          for k in range(len(pid2)):
            for j in range(len(sid2)):
              if sid2[j] < 0:
                if pid2[k] == sid2[j]:
                  match_inrange += 1
                  flag = 'match'
                  break
            if flag == 'nomatch':
              pid_nomatch.append(pid2[k])
          	
          for k in range(len(pid_nomatch)):
            for j in range(len(sid3)):
              if pid_nomatch[k] == sid3[j]:
                bad_phot_match += 1
        
          if len(pid2) > 0:
            spec_frac = (match_inrange + bad_phot_match) / len(pid2)
          else:
            spec_frac = 0
            
          #print(i, len(struct_z), mean)
          #print(match_inrange, bad_phot_match, len(pid2))
          #exit()
          
          #spec_frac = len(sid2) / len(pid2)
          # cut to spectral fraction limits
          if spec_frac >= 0.0:       
            if plot:
              plt.figure(num=i, figsize=(24, 12))
              gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 2]) 
              
              plt.suptitle(field.upper() + ' ID#' + str(i) + '; RA, DEC = [%.4f, %.4f], Spec Frac = %.3f'%(mean_ra, mean_dec, spec_frac))
              
              ax1 = plt.subplot(gs[0])
              edge_detect = np.where(struct_no < 0)
              good_detect = np.where(struct_no > 0)
              plt.errorbar(struct_z[good_detect], struct_flux[good_detect], yerr=struct_flux_err[good_detect], fmt='k.', ms=10)  # Data
              plt.errorbar(struct_z[edge_detect], struct_flux[edge_detect], yerr=struct_flux_err[edge_detect], fmt='r.', ms=10)  # Data
              points = np.linspace(mean-1,mean+1,2000)
            
            posrange = np.where((pfield == field) & (pz < mean + 0.05) & (pz > mean - 0.05))
            
            pras = pra[posrange]
            pdecs = pdec[posrange]
            pzs = pz[posrange]
            name = pname[posrange]
         
            flag = 'new'
            if len(pzs) > 0:
              d = np.empty(len(pzs))
           
              for k in range(len(pzs)):
                d[k] = math.acos(math.sin(pdecs[k]) * math.sin(np.deg2rad(mean_dec)) + math.cos(pdecs[k]) * math.cos(np.deg2rad(mean_dec)) * math.cos(pras[k]-np.deg2rad(mean_ra)))
              d = d * dA
              closest = np.where(d == np.min(d))[0][0]
              
              #print(closest, d[closest], name[closest])
              if d[closest] < 3.0: # Mpc matching radius
                flag = str(name[closest])
                if plot:
                  plt.plot(points, gaus(points, *popt), 'g--', linewidth=3, label='Fit')
                  plt.title('Struct ' + flag + ', Dist: %.2f Mpc'%(d[closest]))
                flag = 's' + flag
                candidate = False
            
            if plot:
              if flag == 'new':
                plt.plot(points, gaus(points, *popt), color='orange', linestyle='dashed', linewidth=3, label='Fit')
                plt.title('Gauss Amplitude = %.1f'%(gauss_amp))
              plt.ylim(-0.05*np.max(struct_flux),1.1*np.max(struct_flux))
              plt.xlim(np.min(struct_z)-0.05, np.max(struct_z)+0.05)
              plt.xlabel('z')
              plt.ylabel('Iso Flux')
              plt.legend(['z = %.4f$\pm$%.4f\n$\sigma_{z}$ = %.4f$\pm$%.4f'%(mean, mean_err, std, std_err)], loc='upper right', fontsize=16)
            
              # histogram of galaxies in spatial distance range
              hist_z = spec_z1
              hist_z_phot = phot_z1
              
              plt.subplot(gs[1])
              y, y2 = 0, 0
              histrange = [mean-delta_redshift,mean+delta_redshift]
              
              y2, bins, patches2 = plt.hist(hist_z_phot, 15, facecolor='lightskyblue',edgecolor='black',linewidth=1, alpha=0.75, range=histrange, label='z$_{phot}$')
              # hatched pattern for second histogram
              plt.rcParams['hatch.linewidth'] = 2.0
              y, bins, patches = plt.hist(hist_z, 15, facecolor='none',edgecolor='#045c5a',linewidth=3,alpha=1, range=histrange, hatch='/', label='z$_{spec}$')
              plt.xlabel('z')
              plt.ylabel('Count')
              plt.title('Objects within 2 Mpc')
              plt.axvline(x=mean-3*std,color='orange',linestyle='dashed', linewidth=3, label='3$\sigma$ Boundary')
              plt.axvline(x=mean+3*std,color='orange',linestyle='dashed', linewidth=3)
              plt.xlim(mean-delta_redshift,mean+delta_redshift)
              peak = np.max([np.max(y2),np.max(y)])
              plt.ylim(0,int(1.2*peak)+1)
              plt.legend(loc='upper right')
              
              # plot galaxies in redshift+distance range
              ax3 = plt.subplot(gs[2])
              ss = np.where((spec_z > mean - delta_redshift) & (spec_z < mean + delta_redshift))
              sp = np.where((phot_z > mean - delta_redshift) & (phot_z < mean + delta_redshift))
              plt.scatter(np.rad2deg(phot_ra[sp]),np.rad2deg(phot_dec[sp]),s=110,edgecolors='k',color='lightskyblue',label='z$_{phot}$')
              plt.scatter(np.rad2deg(spec_ra[ss]),np.rad2deg(spec_dec[ss]),s=120,linewidth=2,edgecolors='#045c5a',facecolor='none',marker='d',label='z$_{spec}$')
              plt.scatter(mean_ra,mean_dec,s=200,edgecolors='k',color='orange',marker='*',label='Centroid')
              # plot circle of 2 Mpc = convert degrees to angular distance
              circ = plt.Circle((mean_ra,mean_dec), 2*radius_Mpc, color='orange', fill=False, linestyle='dashed', linewidth=3, label='2 Mpc Radius')
              ax3.add_artist(circ)
              ax3.ticklabel_format(useOffset=False)
              
              nticks = radius_Mpc*6.0/4.0
              plt.xlim(mean_ra+4.3*radius_Mpc,mean_ra-4.3*radius_Mpc)
              plt.ylim(mean_dec-4.3*radius_Mpc,mean_dec+4.3*radius_Mpc)
              ax3.xaxis.set_ticks(np.arange(mean_ra-4.3*radius_Mpc,mean_ra+4.3*radius_Mpc, nticks))
              ax3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
              plt.title('Objects $\Delta$z < %.2f'%(delta_redshift))
              plt.xlabel('RA ($\degree$)')
              plt.ylabel('Dec ($\degree$)')
              plt.legend(loc='upper right', scatterpoints = 1)
              
              outfile = directory + field + '.' + str(i) + '.png'
              
              plt.subplots_adjust(top=0.90)
              plt.tight_layout(w_pad=0.4)
              plt.savefig(outfile)
              plt.close('all')
        
            print('db01', lk, field, 'd' + str(detect) + '_a' + marea, i, len(struct_z), mean, mean_err, std, std_err, mean_ra, mean_dec, errs, gauss_area, gauss_amp, flag, spec_frac, this_edge_frac)
