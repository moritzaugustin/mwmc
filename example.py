'''
Created on 05.03.2013

@author: moritz
'''

import argparse
import os
import shutil
import tables
import numpy
import subprocess
import socket

# parameters
numtype = 'float' # int or float
nspecies = 100 #179
nrows =  100 #2760
ncols = 150 #4800

# host related stuff
if socket.gethostname() == 'mali':
    python_interpreter = '/opt/anaconda/bin/python2'
    test_folder = '/tmp/mwmc/exampledata'
else:
    python_interpreter = 'python2'
    test_folder = '/lustre/augustin/mwmc/exampledata'

def main():
    parser = argparse.ArgumentParser(description='Generates a couple of test datasets.')
    parser.add_argument('--folder', default='/tmp/mwmc/exampledata', 
                        help='target folder of the testset (default: %(default)s); will be overwritten')
    args = parser.parse_args()
    
    data_folder = test_folder+'/data'
    result_folder = test_folder+'/result'
    
#     print('folder '+test_folder+' will be overwritten. continue? y or n. def: n')
#     if raw_input()!='y':
#         print('exiting.')
#         exit()
    
    try:
        shutil.rmtree(test_folder, ignore_errors=True)
        os.makedirs(data_folder)
        os.makedirs(result_folder)
    except:
        pass
    
    NODATA_value = -9999
    yllcorner = 0.0 # will be ignored
    xllcorner = 0.0 # will be ignored
    cellsize = 0.008333333 # will be ignored
    
    print('generating random test data...')
    for k in range(1, nspecies+1):
        tfile = tables.openFile(data_folder+'/'+'species'+str(k)+'.h5', 
                                mode='w')
        tfile.createArray(tfile.root, 'nrows', nrows)
        tfile.createArray(tfile.root, 'ncols', ncols)
        tfile.createArray(tfile.root, 'NODATA_value', NODATA_value)
        tfile.createArray(tfile.root, 'yllcorner', yllcorner)
        tfile.createArray(tfile.root, 'xllcorner', xllcorner)
        tfile.createArray(tfile.root, 'cellsize', cellsize)
        if numtype == 'int':
            data = numpy.random.randint(0, 2, size=(nrows, ncols)).astype(numpy.int8)
        elif numtype == 'float':
            data = numpy.random.rand(nrows, ncols).astype(numpy.float64)*0.1
        else:
            print('error: numtype {n} unknown. exiting'.format(n=numtype))
            exit()
            
        # some nodata values in each row (and each species)
        for r in range(nrows):
            nodata_inds_col = numpy.unique(numpy.random.randint(0, ncols, 
                                   size=numpy.random.randint(1, ncols)))
            data[r, nodata_inds_col] = NODATA_value
        
        tfile.createArray(tfile.root, 'data', data)
        tfile.close()
        
    print('computing example measure...')
    if numtype=='float':
        measure = 'alpha_shannon_center_unormalized' # alpha_shannon_center, alpha_simpson_center
    elif numtype=='int':
        measure = 'alpha_count_center'
    subprocess.call(python_interpreter+' mwmc.py'+
                    ' --compute'+
                    ' --binfolder='+data_folder+
                    ' --targetfolder='+result_folder+
                    ' --numberformat='+numtype+
                    ' --measure='+measure+
                    ' --processes=3', 
                    shell=True)
    
    print('visualizing result...')
    subprocess.call(python_interpreter+' mwmc.py'+
                    ' --visualize'+
                    ' --targetfolder='+result_folder+
                    ' --numberformat='+numtype,  
                    shell=True)
    
    print('done.')
    

if __name__=='__main__':
    main()
