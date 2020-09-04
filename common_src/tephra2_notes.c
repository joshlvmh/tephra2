
void tephra_calc(ERUPTION *erupt, POINT *pt, WIND *day, STATS *stats, GRAIN *gr) { /* tephra_calc starts ... */
 
   /**********************************************************************************
   * WIND structure:
   * day->windspeed: wind speed in m/s
   * day->wind_direction: wind direction in +/- degrees from north
   * day->wind_height: meters above sea level
   *
   * See common_structures_strat.h for structure definitions.
   **********************************************************************************/
  
  
  int i, j=0; /*, bin = -1; */
  double new_xspace, new_yspace, cos_wind = 0.0, sin_wind = 0.0, windspeed = 0.0;
  double sigma, demon2=0.0, demon3=0.0, ash_fall = 0.0, layer, fall_time_adj = 0.0, total_fall_time=0.0;
  double average_windspeed_x, average_windspeed_y, average_wind_direction, average_windspeed =0.0;
  double x_adj = 0.0, y_adj = 0.0;
  static double min=10e6, max=0.0;
  
#ifdef _PRINT
  fprintf(log_file, "IN tephra_calc ...");
   fflush(log_file);  
#endif

  /* Initialize mass to zero */
  pt->calculated_mass = 0.0;
	
  /* Transform the volcano location coordinate to 0,0 */
  new_xspace = pt->northing - erupt->vent_northing;
  new_yspace = pt->easting - erupt->vent_easting;
/* fprintf(log_file, "point[x][y]=%.0f, %.0f\n",new_xspace, new_yspace); */

  /* do the double integration over grainsize and column height */
  
  /* Interpolate to fine the wind speed and direction below the height of the vent.
     Find one average speed and direction between vent and grid elevation point. 
   * The first values in the wind array give the wind speed and direction at the
   * vent height.
   */
  layer = erupt->vent_elevation - pt->elevation;
  windspeed = (day[0].windspeed * pt->elevation) / erupt->vent_elevation;
  cos_wind = cos(day[0].wind_direction) * windspeed;
	sin_wind = sin(day[0].wind_direction) * windspeed;
  //__assume_aligned(pt, 64);
  //__assume_aligned(pt->calculated_phi, 64);
  #pragma vector aligned
  for (i = 0; i < PART_STEPS; i++) { /* PART_STEPS_LOOP */
    fall_time_adj = 0.0;
    /* Accumulate the particle sizes into bins of whole numbered phi sizes 
    if (!(i % 10)) {
      bin++;
	Initiialize new phi accumulator to zero */
	pt->calculated_phi[i] = 0.0;
#ifdef _PRINT
    fprintf(log_file, "PART_STEP=%d phi[%d] = %g\n", i, i, pt->calculated_phi[i]);
     fflush(log_file);  
#endif
   /* } */
  
    /* Adjust the total fall time of each particle size (i) 
       by the time it takes to descend from vent height to the grid cell (pt) 
     */           
     if (layer > 0) {
     	 fall_time_adj = 
     	 part_fall_time(erupt->vent_elevation, layer, T[i][0].ashdiam, T[i][0].part_density);
     
    }

     for (j = 0; j < COL_STEPS; j++) { /* COL_STEPS_LOOP */
     
    	total_fall_time = T[i][j].total_fall_time + fall_time_adj;
	   
	   /* Calc. x-component and y-component windspeed adjustments  
    	 * for each particle size  falling from each level.
    	 */
     /* removed 2 lines, 2-22-2011
   		x_adj = cos_wind * fall_time_adj * windspeed;		
	  	y_adj = sin_wind * fall_time_adj * windspeed;
	  	*/
	  	/* change 2 lines, 2-22-2011 */
	  	x_adj = cos_wind * fall_time_adj;
	  	y_adj = sin_wind * fall_time_adj;
	  	
			/* Now add the adjustments to the already summed
	     * wind components and
    	 * find the average windspeed in the x and y directions 
    	 * over the total fall time for each particle size at each level.
			*/
			average_windspeed_x = 
			(T[i][j].wind_sum_x + x_adj)/total_fall_time;
			
			average_windspeed_y = 
			(T[i][j].wind_sum_y + y_adj)/total_fall_time;
   	    	
    	/* If zero, make windspeed a very small value (cannot divide by zero in next step) */
      if (!average_windspeed_x) average_windspeed_x = .001;
      if (!average_windspeed_y) average_windspeed_y = .001;
      
      /* Find the average wind direction (direction of the velocity vector) */
    	if (average_windspeed_x < 0) {
    	  average_wind_direction =
    	  atan(average_windspeed_y/average_windspeed_x ) + pi;
    	}
    	else 
    	  average_wind_direction = atan(average_windspeed_y/average_windspeed_x);
    	
    	/* Find the average wind speed ( magnitude of the velocity vector) */
      average_windspeed =
      sqrt(average_windspeed_x * average_windspeed_x + 
      average_windspeed_y * average_windspeed_y);

    	if (total_fall_time > max) max = total_fall_time;
    	if (total_fall_time < min) min = total_fall_time;
    	
    	/* calculate the value of sigma (dispersion) based on total_fall_time  
    	    to acct for the change in the shape of the column with ht - increasing radius
          ht_above_vent = T[i][j].particle_ht - erupt->vent_elevation; 
       */
       
      /* falltime for fine particles */
      if (total_fall_time >= erupt->fall_time_threshold) {      
         sigma = 
         erupt->eddy_constant * pow((total_fall_time + T[i][j].plume_diffusion_fine_particle), 2.5);
      } 
      
      /* falltime for coarse particles */  
      else {
         sigma =
         4.0 * erupt->diffusion_coefficient * (total_fall_time + T[i][j].plume_diffusion_coarse_particle);
      }     
    
      demon2 =  pi * sigma;

      /* Modify fall time by the variation of wind velocity with height */
      demon3 = strat_average( average_wind_direction, 
                       average_windspeed,             
			                 new_xspace, new_yspace, 
			                 total_fall_time,
			                 sigma); 
			         
			 ash_fall = (T[i][j].demon1 / demon2) * demon3;
			 pt->calculated_mass += ash_fall;
			 pt->calculated_phi[i] += ash_fall;
		}  /* COL_STEPS_LOOP */  
   #ifdef _PRINT   
fprintf(log_file, "bin[%g] mass[%g] part[%g]", pt->calculated_phi[i], pt->calculated_mass, ash_fall); 
			fprintf(log_file, "\n"); 
 fflush(log_file);  			 

  #endif
