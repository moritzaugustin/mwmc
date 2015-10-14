'''
@author: Moritz Augustin
'''

import argparse
import subprocess

data_tuples = [('d1_simple', 'float', 'alpha_shannon_center_percellnormalized'),
               ('d1_simple', 'float', 'alpha_shannon_center_unnormalized'),
               ('d2_improved', 'float', 'alpha_shannon_center_percellnormalized'),
               ('d2_improved', 'float', 'alpha_shannon_center_unnormalized'),
               ('d3_bin', 'int', 'alpha_count_center')]

subsets = [None, # no subset = all species
           'subset_ed.txt', 
           'subset_ed_ug.txt', 
           'subset_ed_wa.txt', 
           'subset_edge.txt',
           'subset_edge_ug.txt',
           'subset_edge_wa.txt',
           'subset_rl.txt',
           'subset_rl_ug.txt',
           'subset_rl_wa.txt',
           'subset_ug.txt',
           'subset_wa.txt'
           ]

data_path = '/mnt/storage/augustin/Biodiversity_data'
subsets_path = '/mnt/storage/augustin/Biodiversity_data/subsets_unweighted'
results_path = '/mnt/storage/augustin/Biodiversity_data/results_alpha_center_20151013'
python_interpreter = 'python2' #'/opt/anaconda/bin/python2'
no_procs = 6 # no of parallel processes
no_blocks = 12



# you can run the experiment by specying the data tuple id (i.e. 0, 1, 2, ...)
# via command line => no automatic iteration over alle tuples but 5 seperate jobs possible
# default: iterate over all data_tuples
parser = argparse.ArgumentParser()
parser.add_argument('--data_tuple_id', #required=True, 
                    help='The data tuple id (i.e. 0, 1, ..., '+str(len(data_tuples)-1)+')')
args = parser.parse_args()
data_tuple_id = int(args.data_tuple_id) if args.data_tuple_id else None

all_data_tuple_ids = range(len(data_tuples))

if data_tuple_id not in all_data_tuple_ids:
    data_tuple_ids = all_data_tuple_ids
    print('iterating over all data tuples as --data_tuple_id is not specified')
else:
    data_tuple_ids = [data_tuple_id]

for dtid in data_tuple_ids:
    
    
    print(('experiment started -- data tuple {datatuple} selected (of {nodatatuples} in total)'+
           ' and {nosubsets} subsets'+
           '...\n').format(nodatatuples=len(data_tuples), 
                           nosubsets=len(subsets),
                           datatuple=dtid))
    
    dataset, dataformat, measurename = data_tuples[dtid]
    print('dataset: '+dataset)
    print('dataformat: '+dataformat)
    print('measurename: '+measurename+'\n')
    
    for subset in subsets:
        if subset: 
            subsetparam = ' --subsetfile='+subsets_path+'/'+subset
            subsetresult = '/subset_'+subset
        else: # None -- corresponds to no subset (but all species)
            subsetparam = ''
            subsetresult = ''
    
        command = (python_interpreter+' mwmc.py'+
                        ' --compute'+
                        ' --binfolder='+data_path+'/'+dataset+'/hdf5'+
                        subsetparam+
                        ' --targetfolder='+results_path+'/'+dataset+'/'+measurename+subsetresult+
                        ' --numberformat='+dataformat+
                        ' --measure='+measurename+
                        ' --processes='+str(no_procs)+
                        ' --no_blocks='+str(no_blocks))
        print(command)
    
        subprocess.call(command, shell=True)


print('experiment done.')