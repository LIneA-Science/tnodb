#Filter_Freq_Analysis v1.0
# This program filters the file generated by the Freq_AnaLysis program with options of magnitude
#position error, number of days between exposition, and if there is more than one observation
#in the same night with different filters.
# As output you can generate a table with this data, or just a file with the name of the objects,
#to enter the program Download_from_DES (not yet available)
#
#DAT File
#The dat file must be Like: And be called Filter_Freq_Analysis.dat
'''
Value					| Row name		 	 |Description 
InputFile.txt			| InputName			 |Input Name			
0					| FClass	 	 	 |Class (0 if All DES File, else readme)
55.0					| FVmag	 	 	 |Visual Mag(0 if no filter)(dot decimal separation)
0					| Fposerr 	 	 	 |Position error (0 if no filter)
0					| FObjectName 	 	 |Object name ( 0 if all objects)
0					| FNumNights		 |Minimun Observed nights
0					| ShowMTFS			 |Show Only if Has more than one filter in the some night
0					| FilterDelta		 |Filter my minimun Delta time in days
0					| OnlyNames			 |Output files with only names
OutputFile.txt			| OutputName 		 |Output File name 
'''

#import Libraries
print('importing Libraries')
from astropy.table import Table
from numpy import *

# read Filters file 
print('Reading filter files')
filters = Table.read("Filter_Freq_Analysis.dat",format='ascii', delimiter='|' )

#Read Filters
InputName = filters['Value'][0]
FClass= filters['Value'][1]       
FVmag= filters['Value'][2]
FVmag = float(FVmag)
FPoserr= filters['Value'][3]    
FPoserr=float(FPoserr)
FObjectName= filters['Value'][4]
FNumNights=filters['Value'][5]
FNumNights=int(FNumNights)
ShowMTFS=filters['Value'][6]
FilterDelta =filters['Value'][7]
FilterDelta = int(FilterDelta)
OnlyNames = filters['Value'][8]
OutputName= filters['Value'][9]

#Read DES_FREQ_FILE
print('Reading', InputName, 'files')
freq_des = Table.read(InputName,format='ascii', delimiter='\t' )
freq_des

#apply Filters
if FClass != '0':
    freq_des = freq_des.group_by ("Class")
    mask = freq_des.groups.keys['Class'] == FClass
    freq_des  = freq_des.groups[mask]

if FVmag != 0:
    magnitude = freq_des["VMagMax"] <= FVmag
    freq_des = freq_des[magnitude]
    
if FPoserr != 0:
    position = freq_des["MaxPosErr"] <= FPoserr
    freq_des = freq_des[position]
    
if FObjectName != '0':
    freq_des = freq_des.group_by ("ObjectName")
    mask = freq_des.groups.keys['ObjectName'] == FObjectName
    freq_des  = freq_des.groups[mask]

if FNumNights != 0:
    diffnights = freq_des['DiffNights'] >= FilterDelta
    freq_des = freq_des[diffnights]
    
if ShowMTFS == '1':
    ShowMTFS = 'Y'
    freq_des = freq_des.group_by("MfilterN")
    mask = freq_des.groups.keys['MfilterN'] == ShowMTFS
    freq_des  = freq_des.groups[mask]
    
if FilterDelta != 0:
    deltatime = freq_des['DeltaTime'] >= FilterDelta
    freq_des = freq_des[deltatime]
#Save File
print('Saving', OutputName)
if OnlyNames == '0':
    freq_des = freq_des.group_by("ObjectName")
    freq_des = freq_des.groups.keys
    freq_des.write(OutputName, format='ascii', delimiter='\t')
else:
    freq_des.write(OutputName, format='ascii', delimiter='\t')
#FINISHED
print('finished')
