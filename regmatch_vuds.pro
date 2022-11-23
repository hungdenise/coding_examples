function cosmo, z
	
	omega_m = 0.27
	omega_lambda = 0.73
	return, 1.0/sqrt(omega_m*(1+z)^3 + omega_lambda)

end

pro regmatch_vuds, max_distance
; Match objects detected by SExtractor across different redshift slices based on user set physical separation (Mpc)
; max_distance = max separation between objects for matching in Mpc
	close, /all	

	H0 = 70.0 ; Hubble constant
	sigmas = ['3','3.5'] ; DETECT_THRESH sigma in SExtractor
	marea = '100' ; DETECT_MINAREA in SExtractor
	dbs = ['_db01_', '_db03_'] ; DEBLEND_MINCONT in SExtractor (0.1 or 0.3)
	fields = ['CFHTLS-D1', 'COSMOS', 'ECDFS']

	for f = 0, 2 do begin 
	field = fields[f]
	
	bins = 832
	fits = 'sav_maps/' + field + '_overdens80.fits'

	for ss = 0, 1 do begin
	sigma = sigmas[ss]
	for dd = 0, 1 do begin
	db = dbs[dd]
	
	empty = make_array(bins, /integer)
	bindex = make_array(bins*2, /integer)
	
	outfile = 'regmatch/regmatch_' + field + '_d' + sigma + '_a' + marea + db + strcompress(string(max_distance),/remove_all) + 'd.txt'
	openw, lun1, outfile, /get_lun, width=240
	printf, lun1, 'bin,no,z,area,flux,flux_err,ra,dec,d'
	
	a = 1
	for s = 0, bins-1 do begin
	
	regFile = 'rms/' + field + '_d' + sigma + '_a' + marea + db + strcompress(string(s),/remove_all) + '.cat' ; SExtractor detections, listed in individual redshift slice files 
	if file_lines(regFile) gt 15 then begin 
		
		number2 = 0
		flux2 = 0
		fluxerr2 = 0
		ra2 = 0
		dec2 = 0
		
		readcol,regFile,number2,flux2,fluxerr2,background2,threshold2,area2,ra2,dec2, format='I,D,D,X,X,D,D,I,X,X,D,D,X,X,X',/silent

		if a eq 1 then begin
			number = number2
			flux = flux2
			fluxerr = fluxerr2
			area = area2
			ra = ra2
			dec = dec2
		endif else begin
			number = [number, number2]
			flux = [flux, flux2]
			fluxerr = [fluxerr, fluxerr2]
			area = [area, area2]
			ra = [ra, ra2]
			dec = [dec, dec2]
		endelse
		bindex[a] = bindex[a-1] + n_elements(number2)
		a = a + 1
	endif else begin
		empty[s] = 1
	endelse
	endfor
	
	z1 = make_array(bins)
	z2 = make_array(bins)
	z0 = make_array(bins)
	z = make_array(bins,3)
	
	struct_n = -1
	
    for n = 0, bins-1 do begin  
		h = headfits(fits,ext=n)
		
		z1[n] = sxpar(h,'z1')
		z2[n] = sxpar(h,'z2')
		
		z0[n] = (z1[n]+z2[n])/2.0		; z value of each individual bin
		
		if (n ne 0) then begin
			z[n,0] = (z1[n]+z2[n-1])/2.0		; average z from end of last bin and start of next bin
		
			DC = 299792./H0*QROMB('cosmo', 0.0, z[n,0])
			z[n,1] = DC
		
			DA = z[n,1]/(1+z[n,0])
			z[n,2] = DA
		endif
	endfor
	
	b = -1
	
	for n = 0, bins-2 do begin
		if (empty[n] eq 0) then begin
		
			b = b + 1
			for j = bindex[b], bindex[b+1]-1 do begin
				ra1 = ra[j]*!pi/180.
				dec1 = dec[j]*!pi/180.
			
				if (n lt bins-2) and (empty[n+1] ne 1) then begin
					for k = bindex[b+1],bindex[b+2]-1 do begin
						ra2 = ra[k]*!pi/180.
						dec2 = dec[k]*!pi/180.
		
						cosA = sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2)
						A = acos(cosA)		; angular separation in radians
		
						d = z[n+1,2]*A		; angular diameter distance * angular separation

						x = 1
						if (d le max_distance) and (empty[n+2] ne 1) then begin

							struct_n = struct_n + 1
							
							line1 = strcompress(string(n+1)) + ',' + strcompress(string(number[j])) + ',' + strcompress(string(z0[n])) + ',' + strcompress(string(area[j])) + ',' + strcompress(string(flux[j])) + ',' + strcompress(string(fluxerr[j])) + ',' + strcompress(string(ra[j])) + ',' + strcompress(string(dec[j]))
							line2 = strcompress(string(n+2)) + ',' + strcompress(string(number[k])) + ',' + strcompress(string(z0[n+1])) + ',' + strcompress(string(area[k])) + ',' + strcompress(string(flux[k])) + ',' + strcompress(string(fluxerr[k])) + ',' + strcompress(string(ra[k])) + ',' + strcompress(string(dec[k])) + ',' + strcompress(string(d))
							
							printf, lun1, 'Struct,', struct_n
							printf, lun1, line1
							printf, lun1, line2
						
							d_check = 1
							; continue loop if d < max_distance is found
						
							ra1a = ra1
							ra2a = ra2
							dec1a = dec1
							dec2a = dec2
							flux1a = flux[j]
							flux2a = flux[k]
							flux_tot = flux[j]+flux[k]
							skip = 0

							while (d_check ne 2 or skip eq 0) do begin
								d_check = 2
								
								if (n+x+1 eq bins) then break
								if (empty[n+x+1] eq 1) then break
								for l = bindex[b+x+1],bindex[b+x+2]-1 do begin
									
									; get flux weighted average for linked chain
									ra1 = (flux2a*ra2a+flux1a*ra1a)/flux_tot
									dec1 = (flux2a*dec2a+flux1a*dec1a)/flux_tot
									; replace ra and dec with average of last computed

									ra2 = ra[l]*!pi/180.
									dec2 = dec[l]*!pi/180.
		
									cosA = sin(dec1)*sin(dec2) + cos(dec1)*cos(dec2)*cos(ra1-ra2)
									A = acos(cosA)			; angular separation in radians
			
									d = z[n+x+1,2]*A		; angular diameter distance * angular separation
								
									if (d le max_distance) then begin
										; keep object in bin with smallest distance
										if skip eq 1 then print, 'Skip'
										if d le d_check then begin
											numl = number[l]
											areal = area[l]
											fluxl = flux[l]
											fluxl_err = fluxerr[l]
											ra1a = ra1
											ra2a = ra2
											dec1a = dec1
											dec2a = dec2
											flux1a = flux1a + flux2a
											flux2a = flux[l]
											flux_tot = flux1a + flux2a
											; keep relevant coordinates to match with next iteration on while loop
										endif
										d_check = d
									endif
								endfor
								
								if d_check eq 2 then skip = 1
								
								if d_check ne 2 then printf, lun1, strcompress(string(n+x+2)) + ',' + strcompress(string(numl)) + ',' + strcompress(string(z0[n+x+1])) + ',' + strcompress(string(areal)) + ',' + strcompress(string(fluxl)) + ',' + strcompress(string(fluxl_err)) + ',' + strcompress(string(ra2a*180./!pi)) + ',' + strcompress(string(dec2a*180./!pi)) + ',' + strcompress(string(d_check))
								x = x + 1
							endwhile
						endif
					endfor
				endif
			endfor
		endif
	endfor
	
	free_lun, lun1
	endfor
	endfor
	endfor
end