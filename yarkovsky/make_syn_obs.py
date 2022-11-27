# get ephemerides files from OrbFit output and create sets of synthetic observations of varying arc lengths
import os
import numpy as np

# add noise to an ephemeris position to create a synthetic observation
def add_noise(ephem_ra, ephem_dec, noise_level):
  ra = ephem_ra.split()
  dec = ephem_dec.split()
  
  # convert sexagesimal to degrees; declination treated only as positive for now
  ra_deg = float(ra[0]) * 15 + float(ra[1]) / 60 * 15 + float(ra[2]) / 3600 * 15
  dec_deg = float(dec[0][1:]) + float(dec[1]) / 60 + float(dec[2]) / 3600

  # add Gaussian noise to the ra and dec independently; ra noise is modified by the cos(Dec) term
  new_dec = dec_deg + np.random.normal(0, noise_level, 1)[0] / 3600
  new_ra = ra_deg + np.random.normal(0, noise_level, 1)[0] / 3600 / np.cos(np.deg2rad(new_dec))
  
  # add noise from the precision limits on ra/dec in quadrature with the noise level to get expected uncertainty
  ra_rms_noise = "{0:.3f}".format(np.sqrt(noise_level**2 + (0.0005 * 15 * np.cos(np.deg2rad(new_dec)))**2))
  dec_rms_noise = "{0:.3f}".format(np.sqrt(noise_level**2 + 0.005**2))
  
  # convert degrees back to sexagesimal
  noise_ra_h = int(new_ra / 15)
  noise_ra_m = int(new_ra / 15 % 1 * 60)
  noise_ra_s = "{0:.3f}".format(new_ra / 15 % 1 * 60 % 1 * 60)
  noise_ra_all = [noise_ra_h, noise_ra_m, noise_ra_s]
  
  noise_dec_d = int(new_dec)
  noise_dec_m = int(new_dec % 1 * 60)
  noise_dec_s = "{0:.2f}".format(new_dec % 1 * 60 % 1 * 60)
  noise_dec_all = [noise_dec_d, noise_dec_m, noise_dec_s]
  
  # put together correctly formatted strings
  noise_ra = ''
  noise_dec = ''
  
  for e in noise_ra_all:
    if float(e) < 10:
      noise_ra += '0'
    noise_ra += str(e) + ' '
  
  for e in noise_dec_all:
    if float(e) < 10:
      noise_dec += '0'
    noise_dec += str(e) + ' '
  
  # now put back the original sign on dec
  noise_ra = noise_ra[:-1]
  noise_dec = dec[0][0] + noise_dec[:-1]
  
  return noise_ra, noise_dec, ra_rms_noise, dec_rms_noise

# first generate files to find ephemerides, run OrbFit, then make modified observations file
def make_inputs_for_synobs():
  import os
  from astropy.time import Time
  
  # get the reference OrbFit options file 
  f = open('963.oop', 'r')
  ref_oop = f.read()
  f.close()

  asts = ['45898', '50219', '13936'] # asteroids by designation number
  
  directory = 'inputs/get_prephem/'
  apriori_types = ['_A05', '_A1', '_A2', '_A4', '_A8'] # choice of a priori A2 values

  oep_path = 'inputs/get_ephem/oep/' # OrbFit writes out ephemerides to .oep files
  oeps = os.listdir(oep_path)

  for ast in asts:
    if len(oeps) == 0: # if no ephemerides files, then create options files for generating the ephemerides
      oop = directory + ast + '_apriori' # file prefix for new options file
      new_oop = ref_oop.replace('963', ast)
      
      # turn on .yarko_apriori and set .yarko_std to a very small value
      new_oop = new_oop.replace('.yarko_apriori=.FALSE.', '.yarko_apriori=.TRUE.')
      new_oop = new_oop.replace('.yarko_std=3.25d-4', '.yarko_std=1.0d-25') 
      
      # set the .yarko_nominal to the desired value and write out to the individual options files
      new_oop2 = new_oop.replace('.yarko_nominal=0.d0', '.yarko_nominal=8d-5')
      f2 = open(oop + '_A8.oop', 'w')
      f2.write(new_oop2)
      f2.close()
            
      new_oop2 = new_oop.replace('.yarko_nominal=0.d0', '.yarko_nominal=4d-5')
      f2 = open(oop + '_A4.oop', 'w')
      f2.write(new_oop2)
      f2.close()
      
      new_oop2 = new_oop.replace('.yarko_nominal=0.d0', '.yarko_nominal=2d-5')
      f2 = open(oop + '_A2.oop', 'w')
      f2.write(new_oop2)
      f2.close()
      
      new_oop2 = new_oop.replace('.yarko_nominal=0.d0', '.yarko_nominal=1d-5')
      f2 = open(oop + '_A1.oop', 'w')
      f2.write(new_oop2)
      f2.close()
      
      new_oop2 = new_oop.replace('.yarko_nominal=0.d0', '.yarko_nominal=5d-6')
      f2 = open(oop + '_A05.oop', 'w')
      f2.write(new_oop2)
      f2.close()
      
    else:
      noise_as = [0.03] # noise in arc seconds
      arcs = np.arange(200, 0, -25) # get arc lengths from 25 to 200 years, in 25 year steps
      inityear = 1972 # first year in the ephemerides
      
      # get original observations file
      f = open('inputs/obs/rwo_gdr3_t12/' + ast + '.rwo', 'r')
      rwo_lines = f.readlines()
      f.close()
      
      for na in noise_as:
        print(ast, na)
        for apriori_type in apriori_types:
          # do 100 iterations for each synthetic observations set
          for i in range(100):
            for arc in arcs:
              f_arc = open('inputs/obs/get_ephem_span/' + ast + apriori_type + '_as' + str(na) + '_' + str(i) + '_arc' +str(arc) + '.rwo', 'w')
              for line in rwo_lines[:7]:
                f_arc.write(line) # copy the same header info
              f_arc.close()
            
            line = rwo_lines[-1] # use the last line as a template for a synthetic observation line
            
            # get the ephemeris positions
            fe = open(oep_path + ast + '_apriori' + apriori_type + '.oep') 
            elines = fe.readlines()
            fe.close()
            
            for eline in elines[9:]:
              # convert .oep date to .rwo 
              mjd = Time(int(eline[19:25]), format='mjd')
              oep_date = mjd.iso.split(' ')[0].replace('-',' ') + eline[25:32]
              
              # get ephemeris position ra/dec
              ephem_ra = eline[35:47]
              ephem_dec = eline[49:61]
              
              # get new ra/dec with noise added to ephemeris positions
              noise_ra, noise_dec, ra_rms_noise, dec_rms_noise = add_noise(ephem_ra, ephem_dec, na)
              
              # replace bias/residuals with 0.000
              weightra = ra_rms_noise + ' T    0.000    0.000'
              weightdec = dec_rms_noise + ' T    0.000    0.000'
              
              # synthetic observation line for the .rwo file
              new_line = line[:17] + oep_date + line[34:45] + 'E-06 ' + noise_ra + line[62:69] + 'E-02' + line[73:77] + weightra + line[102:103] + noise_dec + line[115:122] + 'E-02' + line[126:130] + weightdec + line[155:178] + 'V 500' + line[183:]
              
              # check the year and only append observation in .rwo file if before arc length cutoff
              year = int(oep_date[:4])
              
              for arc in arcs:
                if year <= inityear + arc:
                  f = open('inputs/obs/get_ephem_span/' + ast + apriori_type + '_as' + str(na) + '_' + str(i) + '_arc'  + str(arc) + '.rwo', 'a')
                  f.write(new_line)
                  f.close()

if __name__=='__main__':
  make_inputs_for_synobs()
