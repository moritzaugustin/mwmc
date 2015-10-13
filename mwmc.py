'''
Initial version written in 2013.
Updated: 2015-10-13

@author: Moritz Augustin
'''

# modules from python std library
import argparse
import os
import re
import time
import multiprocessing
import subprocess

# third party modules
import numpy
import tables
import numba

# own modules
import measures

# GLOBAL PARAMETERS
# computation parameters
window_size_default = 3
no_processes_default = 1

# files, data format parameters
srcfile_suffix = '.asc'
src_folder_default = '../data'
nodata_global = -1 # the value NODATA will be replaced with
appended_cells_value = nodata_global # value, for the grid cells outside the grid (left/right/bottom/top)
binfile_suffix = '.h5'
bin_folder_default = '../bindata'
binfile_numpy_dtype = None # one of the two  binfile types below will be chosen by commandline
binfile_numpy_dtype_float = numpy.float64
binfile_numpy_dtype_int = numpy.int8
resfile_numpy_dtype = numpy.float64
target_folder_default = '../results/'+time.strftime('%Y-%m-%d_%H:%M:%S')
window_size_max = 10003 # (arbitrary chosen) maximal allowed window size 

def convert(args):
    no_conv=0
    print('converting all files '+'*'+srcfile_suffix+' in folder '+args.srcfolder)
    
    # delete all binfiles, such that deleted species will not become ghosts later but are 
    # just gone
    for item in os.listdir(args.binfolder):
        if re.search(binfile_suffix+'$', item):
            os.unlink(args.binfolder+'/'+item)
    
    # convert the files
    for item in os.listdir(args.srcfolder):
# sortiere nach nummer dateinamen... unvollstaendig...
#                   sorted, cmp=lambda x,y: int(re.sub('([0-9])*', r'\1', x))-
#                                           int(re.sub('([0-9])*', r'\1', y))):
#         if re.search('[0-9]+'+srcfile_suffix, item):
        if re.search(srcfile_suffix+'$', item):
            print('  processing '+item)
            bin_filename = args.binfolder+'/'+item.replace(srcfile_suffix, binfile_suffix)
            convert_single(src_filename=args.srcfolder+'/'+item, 
                           bin_filename=bin_filename)
            print('  written into '+bin_filename)
            no_conv += 1
        else:
            print('  skipping '+item)

def convert_single(src_filename, bin_filename):
    """ Converts a file from human-readable ascii format ASC into binary format HDF5.
    
    All values which are NODATA_value will replaced by the value of nodata_global. 
    Within the HDF5 file typical header parameters from the ASC file are saved, in 
    particular: ncols, nrows, NODATA_value, xllcorner, yllcorner, cellsize
    
    Args:
        src_filename (str): The filename of the source asc file.
        bin_filename (str): The filename of the target h5 file.
    """
    
    # define header/data flags and dict for variables defined in asc file
    header = True
    data = False
    variables = {}
    
    # open the file
    fh_src = open(src_filename)
    
    # move through all lines in the file
    for line in fh_src:
        split_list = line.split()
        # we are still in the header
        if header:
            
            # a header row
            if len(split_list) == 2:
                varname, val = split_list
                if varname in ['ncols', 'nrows', 'NODATA_value']:
                    variables[varname] = int(val)
                elif varname in ['xllcorner', 'yllcorner', 'cellsize']:
                    variables[varname] = float(val)
                else:
                    print('Warning: Unknown header variable '+varname+'='+val)
                    variables[varname]
                    
            # this corresponds to the row between header and data
            else:
                data_array = numpy.zeros((variables['nrows'], variables['ncols']), dtype=binfile_numpy_dtype)
                data_rownr = 0 # index in y direction
                header = False
                data = True
                
        # do not use else because otherwise the first row of the data would be lost
        if data:
            # initialize with a value -2, which should be overwritten everywhere
            # use here normal int type to allow for large negative NODATA values
            curr_row_array = -2*numpy.ones(variables['ncols']) 
            
            # iterate over all columns for the current row
            for j, elem in enumerate(split_list):
                # last element contains a newline, strip it!
                if j==len(split_list)-1:
                    elem = elem.rstrip()
                    
                
                # make number fit into data type binfile_numpy_dtype
                # however too large negative values are not checked 
                # and only the NODATA_value case is considered
                if binfile_numpy_dtype==numpy.float64:
                    elem_num = float(elem)
                else:
                    elem_num = int(elem) # convert string to integer
                    dtype_max = numpy.iinfo(binfile_numpy_dtype).max
                    if elem_num>dtype_max:
                        elem_num = dtype_max
                        print('Warning: In file '+src_filename+' are data values which exceed the'+
                              ' maximal possible value ('+str(dtype_max)+') of the current data type')
                        
                curr_row_array[j] = elem_num
                
            # set NODATA_values (can exceed data type) to nodata_global (must not exceed data type)
            curr_row_array[curr_row_array==variables['NODATA_value']] = nodata_global
            data_array[data_rownr,:] = curr_row_array.astype(binfile_numpy_dtype)
            data_rownr += 1
    
    fh_src.close()
    
#     print("converted {src}".format(src=src_filename))
#     print(data_array)
    
    # create and save h5 file
    h5f = tables.openFile(bin_filename, 'w')
    
    # write parameters
    for varname, val in variables.iteritems():
        h5f.createArray('/', varname, val)
        
    # write data
    h5f.createArray('/', 'data', data_array)
    
    h5f.close()

def compute(args):
    
    # todo
    # 1) preparation??
    
    # todo
    # generate HPP file (opt) 
    
    # todo
    # compile (opt)
    
    
    print('starting computation of the measure on the full grid')
    
    # collecting data for any case: multiproc, singleproc
    
    # consider subsets or all species
    if args.subsetfile:
        subsetfile = str(args.subsetfile)
        try:
            f = open(subsetfile, 'r')
            print('subsetfile could be opened: '+subsetfile)
            
            if not args.weighted_subset:
                species_subset  = [line.strip() for line in f]
                binfile_list = [item for item in os.listdir(args.binfolder) if 
                                re.search(binfile_suffix+'$', item) and 
                                item in [spec+binfile_suffix for spec in species_subset]]
            else:
                species_subset  = [line.strip() for line in f]
                binfile_list = []
                for s in species_subset:
                    species, count = s.split()
                    for i in range(int(count)):
                        item = species+binfile_suffix
                        if item in os.listdir(args.binfolder):
                            binfile_list.append(item)
        except IOError: 
            print('error: subsetfile '+subsetfile+' does not exist')
            exit()
        else: 
            f.close()
    else:
        print('kein subsetfile gesetzt')
        binfile_list = [item for item in os.listdir(args.binfolder) if 
                        re.search(binfile_suffix+'$', item)] 
    
    no_species = len(binfile_list)
    if no_species == 0:
        print('error: no species files found - did yoy forget to convert the data?')
        exit()
    
    # read grid parameters from the first file in the list
    h5f = tables.openFile(args.binfolder+'/'+binfile_list[0])
    no_rows = h5f.root.nrows.read()
    no_cols = h5f.root.ncols.read()
    h5f.close()
    
    # set row_start and row_count to point to the full grid if not differently specified 
    if args.rowstart >= 0:
        row_start_compute = args.rowstart
    else:
        row_start_compute = 0
        
    if args.rowcount >= 1:
        row_count_compute = args.rowcount
    else:
        row_count_compute = no_rows
    
    # actual computation branches
    no_procs = min(args.processes, row_count_compute)
    comp_args_list = []
    no_blocks = max(no_procs, args.no_blocks)
    for p in range(no_blocks):
        
        # divide the rows into equal large blocks (only the last block has to be 
        # adapted such as all rows are distributed)
        if p == no_blocks-1:
            row_count_single = row_count_compute - p*(row_count_compute//no_blocks)
        else:
            row_count_single = row_count_compute//no_blocks
        row_start_single = row_start_compute + p*(row_count_compute//no_blocks)
        
        # compute_single  parameter (and return) dictionary
        comp_arg = {}
        comp_arg['row_start'] = row_start_single
        comp_arg['row_count'] = row_count_single
        comp_arg['no_cols'] = no_cols
        comp_arg['no_species'] = no_species
        comp_arg['binfilenames'] = [args.binfolder+'/'+binfile for binfile in binfile_list]
        comp_arg['winsize'] = args.winsize
        comp_arg['measure'] = args.measure
        
        comp_args_list.append(comp_arg)
           
    # matrix where to collect the result rows all processes
    result_matrix = -1*numpy.ones((row_count_compute, no_cols), dtype=resfile_numpy_dtype)  
    
    print(('starting computation with {noprocs} processes').
          format(noprocs=no_procs))
    
    if no_procs > 1:
        # multiproc version
        pool = multiprocessing.Pool(processes=no_procs)
        results_iter = pool.imap_unordered(compute_single, comp_args_list)
    else:
        # single proc version
        results_iter = [compute_single(comp_arg) for comp_arg in comp_args_list]
    
    # collect results (either multi or singleproc
    for l, comp_arg in enumerate(results_iter):
        result_matrix[comp_arg['row_start']:comp_arg['row_start']+comp_arg['row_count'], :] = \
            comp_arg['result']
        print('{l}/{nblocks} blocks finished (took {rtime}; pure comp: {ctime})'.format(
                           l=l+1, 
                           nblocks=no_blocks, # updated 
                           rtime=comp_arg['runtime'],
                           ctime=comp_arg['comptime']))
    
    # close pool of processes in case of multiproc
    if no_procs > 1:
        pool.close()
        
    # write result matrix in files
    # and consider bring this in own function which receives comp_args_lists[0] and result_matrix
    # first h5 then asc
    # todo: do better than this excpetion avoidance trick
    try: 
        os.makedirs(args.targetfolder)
    except:
        pass
    
    h5f = tables.openFile(args.binfolder+'/'+binfile_list[0]) # for some info
    NODATA_value = h5f.root.NODATA_value.read()
    yllcorner = h5f.root.yllcorner.read()
    xllcorner = h5f.root.xllcorner.read()
    nrows = h5f.root.nrows.read()
    ncols = h5f.root.ncols.read()
    cellsize = h5f.root.cellsize.read()   
    h5f.close()
    
    result_h5 = tables.openFile(args.targetfolder+'/result_winsize'+
                                 str(args.winsize)+binfile_suffix, 'w')
    result_h5.createArray('/', 'data', result_matrix)
    result_h5.createArray('/', 'NODATA_value', NODATA_value)
    result_h5.createArray('/', 'yllcorner', yllcorner)
    result_h5.createArray('/', 'xllcorner', xllcorner)
    result_h5.createArray('/', 'nrows', nrows)
    result_h5.createArray('/', 'ncols', ncols)
    result_h5.createArray('/', 'cellsize', cellsize)
    result_h5.close()
    
    
    
    # now asc    
    result_asc = open(args.targetfolder+'/result_winsize'+str(args.winsize)+srcfile_suffix, 'w')
    result_asc.write('ncols {0}\n'.format(ncols))
    result_asc.write('nrows {0}\n'.format(nrows))
    result_asc.write('xllcorner {0}\n'.format(xllcorner))
    result_asc.write('yllcorner {0}\n'.format(yllcorner))
    result_asc.write('cellsize {0}\n'.format(cellsize))
    result_asc.write('NODATA_value {0}\n'.format(NODATA_value))
    
    for i in range(nrows):
        for j in range(ncols):
            if j == ncols-1:
                delim = '\n'
            else:
                delim = ' '
            # write element i,j into file
            result_asc.write('{beta}{delim}'.format(beta=result_matrix[i,j], delim=delim))
    
    result_asc.close()
    
    return result_matrix

def compute_single(comp_arg):
    """ Computes the measure using one process for a contiguous set of rows.
    
    Args:
        comp_arg (dict): contains....
    """
    
    row_start = comp_arg['row_start']
    row_count = comp_arg['row_count']
    no_cols = comp_arg['no_cols']
    no_species = comp_arg['no_species']
    binfilenames = comp_arg['binfilenames']
    winsize = comp_arg['winsize']

    starttime_single = time.time() # time measurement

    # srcdata will contain all data necessary for this block of rows, meaning 
    # top/left/right/bottom are added (winsize-1)/2 cells,
    # which leads to different indices! 
    add_len = (winsize-1)/2
    data_local = numpy.zeros( (row_count+2*add_len, 
                                    no_cols+2*add_len, 
                                    no_species), 
                                  dtype=binfile_numpy_dtype)
    result_local = numpy.zeros((row_count, no_cols), dtype=resfile_numpy_dtype)
    species_id = -1
    for binfilename in comp_arg['binfilenames']:
        h5f = tables.openFile(binfilename)
        data_global = h5f.root.data.read()
        species_id += 1 
        
        # old: 
        #-1 + int(re.search('([0-9]+)'+binfile_suffix+'$', binfilename).
        #         group(1)) # -1 for species number to species index
        # inner data
        data_local[add_len:-add_len, 
                   add_len:-add_len,
                   species_id] = data_global[row_start:row_start+row_count, :]
                 
        # left and right border cells get appended with default value
        data_local[:,:add_len, species_id] = appended_cells_value
        data_local[:,-add_len:, species_id] = appended_cells_value
          
        # top border cells get appended 
        # first initialize all cells with default values but after that
        # compute the real overlap with the upper data block and get the 
        # data from data_global array
        append_rows = appended_cells_value*numpy.ones((add_len, no_cols+2*add_len))
        no_overlap_rows_top = add_len if row_start >= add_len else row_start
        if no_overlap_rows_top>0: # update only if not first processor
            append_rows[-no_overlap_rows_top:, add_len:-add_len] = \
                data_global[row_start-no_overlap_rows_top:row_start, :]
        data_local[:add_len,:,species_id] = append_rows
          
        # same for bottom border cells
        append_rows = appended_cells_value*numpy.ones((add_len, no_cols+2*add_len))
        # todo must be rechecked by example whether it is correct (i am tired)
        no_rows_global = data_global.shape[0]
        index_lastrow = row_start + row_count
        if no_rows_global - index_lastrow > add_len:
            no_overlap_rows_bottom = add_len
        else: 
            no_overlap_rows_bottom = no_rows_global - index_lastrow
        if no_overlap_rows_bottom > 0: # update only if not last processor
            append_rows[:no_overlap_rows_bottom, add_len:-add_len] = \
                data_global[index_lastrow:index_lastrow+no_overlap_rows_bottom]
                # these are the add_len rows after the rows corresponding the cur proc
        data_local[-add_len:,:,species_id] = append_rows
        
        h5f.close()
    
    comptime_single_start = time.time()

    # call the numba jit-compiled function
    comp_core(data_local, result_local, add_len, 
              measures.measures_dict[comp_arg['measure']])
    
    comp_arg['result'] = result_local
    comp_arg['runtime'] = time.time() - starttime_single
    comp_arg['comptime'] = time.time() - comptime_single_start    
    
    return comp_arg

@numba.jit
def comp_core(data_local, result_local, ws_h, measure):
    row_count = result_local.shape[0]
    no_cols = result_local.shape[1]
    for i in range(row_count):
        offset_i = ws_h+i
        for j in range(no_cols):
            offset_j = ws_h+j
            data_win = data_local[offset_i-ws_h:offset_i+ws_h+1,
                                  offset_j-ws_h:offset_j+ws_h+1]
            
            result_local[i, j] = measure(data_win)
# comp_core_jit = numba.jit(numba.double[:,:](numba.double[:,:,:]))(comp_core)


def main():
    
    # measure runtime
    time_global_start=time.time()
    
    # PARSE COMMANDLINE
    parser = argparse.ArgumentParser(description='This program computes user-defined '+
                                     'measures in a grid of cells using a moving window. '+
                                     'The computational work is distributed to a number of '+
                                     'parallel processes where each of those '+
                                     'calculates the measure for a subset of the '+
                                     'independent rows of the full grid.', 
                                     epilog='The program was written by Moritz Augustin in '+
                                            '2013 and 2014.')
    
    parser.add_argument('--convert', action='store_true', help='convert the data from '+
                        srcfile_suffix+' into '+binfile_suffix+' (HDF5) format')
    parser.add_argument('--srcfolder', default=src_folder_default, help='folder of the source '+
                        srcfile_suffix+' files (default: %(default)s). The file names for the '
                        +'different species have to end with "1'+srcfile_suffix+'" to '+ 
                        '"M'+srcfile_suffix+'"')
    parser.add_argument('--binfolder', default=bin_folder_default, help='folder of the source '+
                        binfile_suffix+' files (default: %(default)s)')
    parser.add_argument('--numberformat', required=True,
                        help='data format of the input numbers (either int or float)')
    parser.add_argument('--measure', help='the measure which is applied. choose from: '+
                                          str(measures.measures_dict.keys()))
    parser.add_argument('--targetfolder', default=target_folder_default, 
                        help='folder where the result files are written into '+
                             '(default: %(default)s)')
    parser.add_argument('--subsetfile', default=None, 
                        help='file with a subset of speciess with one species file in each line '+
                        '(without extension) ')
    parser.add_argument('--weighted_subset', action='store_true', help='flag if the subset'+ 
                        ' file consists of a weight of each species (basically duplicates '+
                        ' the species).')
    parser.add_argument('--no_blocks', type=int, default=-1,
                        help='decompose all rows into this number of blocks for reducing parallel '+
                        'memory consumption (default: no of processes)')
    
    
    parser.add_argument('--compute', action='store_true', help='compute the measure')
    parser.add_argument('--visualize', action='store_true', help='visualize the measure')
    parser.add_argument('--winsize', default=window_size_default, type=int, 
                        help='Size of the moving window (one edge; default: %(default)s)')
    
    parser.add_argument('--processes', default=no_processes_default, type=int, 
                        help='Number of processes to be started in parallel '+
                        '(default: %(default)s)')
    
    # the following options are intended for internal use only
    
    # compute the measure on the current host using multiple processes
    parser.add_argument('--singlehost', action='store_true', 
                        help=argparse.SUPPRESS)
    # 'first row this host has to compute (counting from 0)'
    parser.add_argument('--rowstart', default=-1, type=int, 
                        help=argparse.SUPPRESS)
    # number of rows this host has to compute
    parser.add_argument('--rowcount', default=-1, type=int, 
                        help=argparse.SUPPRESS) 

    args = parser.parse_args()
    
    if args.winsize not in (2*i+1 for i in xrange((window_size_max-1)/2)):
        print('error: winsize must be an odd number and below '+str(window_size_max))
        exit()
        
    if args.processes < 1:
        print('error: number of processes (={0}) must be a positive integer'.format(args.processes))
        exit()
       
    if args.no_blocks <= 0:
        args.no_blocks = -1
        
    global binfile_numpy_dtype
    if args.numberformat == 'int':
        binfile_numpy_dtype = binfile_numpy_dtype_int
    elif args.numberformat == 'float':
        binfile_numpy_dtype = binfile_numpy_dtype_float 
    else:
        print('error: number format must be either int or float')
    
    # CONVERT DATA INTO BINARY HDF5 FORMAT
    if args.convert:
        convert(args)        
                
    # COMPUTE THE MEASURE ASSUMING ALREADY CONVERTED DATA AND WRITE HDF5 RESULT FILES
    if args.compute:
        if args.measure not in measures.measures_dict:
            print('error: the measure {m} does not exist'.format(m=args.measure))
            exit()
        compute(args)
    
    # VISUALIZE THE MEASURE ASSUMING ALREADY COMPUTED RESULTS IN THE TARGET FOLDER
#     if args.visualize:
#         visualize(args)
# imshow(x, interpolation='nearest')
# colorbar()
# xticks((0,2,4),('a','b','end'))
       
    time_global = time.time() - time_global_start
    print('execution of '+__file__+' took  {t:.1f}s'.format(t=time_global))
       
        
if __name__=='__main__':
    main()