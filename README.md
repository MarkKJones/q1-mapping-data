SHMS Q1 mapping data
---------------------------
Q1 mapping data taken in June 2015.

1. The report is in file : Analysis of JLAB Q1 field plots at JLAB.pdf

2. Additional report in file: Determining the cold mass positions from field map.doc

3. An apparatus for hold the hall probe was built so that 64 measurements could be made around a circle.

4. The apparatus could be moved to 47 different z positions.


## Spreadsheets and dat files

1. A scan along Z was taken at current of 1228.2A . The scan was taken at angle setting number 56 and at all 47 z positions.
  * Spreadsheet : Spreadsheets/map-a56-1228A.xls
  * Data file : Datfiles/q1-1228A-z-scan-ang-56.dat

2.  A scan along Z was taken at current of 2169.1A . The scan was taken at angle setting number 55 and 7 for half 47 z positions.
  * Spreadsheet : Spreadsheets/map-2169A-z-half.xls
  * Data file : Datfiles/q1-2169A-z-scan-ang-55.dat
  * Data file : Datfiles/q1-2169A-z-scan-ang-7.dat

3.  A scan along Z was taken at current of 2454.65A . The scan was taken at angle setting number 55 and 7 for half 47 z positions.
  * Spreadsheet : Spreadsheets/map-2454A-z-half.xls
  * Data file : Datfiles/q1-2454A-z-scan-ang-55.dat
  * Data file : Datfiles/q1-2454A-z-scan-ang-7.dat



## Plotting codes

1. All stored in plot_q1_routines.C  :   .L plot_q1_routines.C++

2. Read_angle_file()  : this reads in ang_index_versus_pos.dat and fill global vector angle

3. Read_zpos_file(): This reads in file that is a table of the z position number and actual position in mm. The method get_pos(num) is used to get the position given the position number ( the position number is what is usually in the spreasheet). The z position number and actual positions was taken from spreadsheet map-a56-1228A.xls

4. Set_probe()  : This fits the hall probe data so the the Q1 field measurements can be corrected. Needed by field_corr() method.


5. Plot_z0_data( fname) : Reads data from fname . 
   * Assumes that it is data which is the full angular coverage at a given Z position 
   * Fits the data with fit_function() which is a sum over i=0 to 4 of A_i*sin(2*(i+1)x) + B_i*cos(2*(i+1)x)
   * So p0=A2 and fills the global vector Q1_central_field_data

6. Plot_q1_bi( Current, measured field/current, tosca field/current)

7. run_plot_q1_bi()
   * call set_probe
   * call Read_angle_file()
   * Loops over the Q1 field measurements done at z=0 and 4 different currents.
   * Call Plot_zo_data to determine Q1_central_field_data for each current
   * Call plot_q1_bi

8. run_plot_q1_zscan(TString fname)
   * call Set_probe
   * reads data to fill vectors with zpos_num and zpos_val
   * call Plot_q1_zscan(fname) which reads in data , corrects field for probe non-linearity and makes plot.
 
