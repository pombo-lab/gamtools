import pandas as pd
import numpy as np
import os
import itertools
import h5py
from multiprocessing import Pool


class HDF5FileExistsError(Exception):
    pass

class GamFrequencyMatrix(object):
    """A class for abstracting access to the HDF5 store underlying a frequency matrix"""
    
    def __init__(self, hdf5_store):
        
        self.data = self.data_from_store(hdf5_store)
        
    def data_from_store(self, hdf5_store):
        """Get a frequency matrix from the store"""
        
        return hdf5_store["processed_data"]["frequencies"]
    
    def any_empty(self, loc1_start, loc1_stop, loc2_start, loc2_stop):
        
        data_array = self.data[loc1_start:loc1_stop, loc2_start:loc2_stop]
        
        for row in data_array:
            for cell in row:
                if not cell.sum():
                   return True
        
        return False
    
    def get_matrix(self, loc1_start, loc1_stop, loc2_start, loc2_stop):
        
        if self.any_empty(loc1_start, loc1_stop, loc2_start, loc2_stop):
            return None
        
        else:
            print 'Using cache'
            return self.data[loc1_start:loc1_stop, loc2_start:loc2_stop]
        
    def cache_freqs(self, loc1_start, loc1_stop, loc2_start, loc2_stop, freqs):
        
        self.data[loc1_start:loc1_stop, loc2_start:loc2_stop] = freqs
        
    @staticmethod
    def create_freq_matrix(hdf5_store, no_windows):
        """Create a frequency matrix dataset stored in the HDF5 file"""
        
        grp = hdf5_store.create_group("processed_data")
        
        freq_matrix = grp.create_dataset("frequencies", (no_windows,no_windows,2,2), dtype='i')
        
        return freq_matrix
    
    @staticmethod
    def from_no_windows(hdf5_store, no_windows):
        
        matrix = GamFrequencyMatrix.create_freq_matrix(hdf5_store, no_windows)
        
        return GamFrequencyMatrix(hdf5_store)

class GamExperimentalData(object):
    """A class for abstracting access to the original experimental segmentation"""
    
    def __init__(self, hdf5_store, pseudocount = 0):
        
        self.store = hdf5_store
        
        self.data = self.data_from_store()
        
        self.no_windows = len(self.data.columns)
        
        self.pseudocount = pseudocount
        
    def data_from_store(self):
        """Retrieve data from an hdf5 store and use it to recreate the pandas DataFrame"""
        
        df_data = np.array(self.store["experimental_data"]["segmentation"][:])
        df_columns = pd.MultiIndex.from_tuples(map(lambda i: tuple(i), np.array(self.store["experimental_data"]["columns"][:])))
        df_index = np.array(self.store["experimental_data"]["index"][:])

        return pd.DataFrame(data=df_data, index=df_index, columns=df_columns)
        
    @staticmethod
    def from_multibam(multibam_file, hdf5_store):
        """Read in a multibam file and return a GamExperimentalData object"""
        
        experimental_data = GamExperimentalData.read_multibam(multibam_file)
        
        GamExperimentalData.data_to_store(experimental_data, hdf5_store)
        
        return GamExperimentalData(hdf5_store)
        
    @staticmethod    
    def read_multibam(input_file):
        """Read in a multibam file and return a Pandas DataFrame"""
        
        return pd.read_csv(input_file, delim_whitespace=True,index_col=[0,1]).transpose()
    
    @staticmethod
    def data_to_store(data, hdf5_store):
        """Save experimental data to an hdf5 store"""
        
        grp = hdf5_store.create_group("experimental_data")
        
        #print data
        segmentation = np.array(data)
        #print segmentation
        segmentation_hd = grp.create_dataset("segmentation", segmentation.shape, segmentation.dtype)
        
        segmentation_hd.write_direct(segmentation)
        
        index = np.array(list(data.index))
        
        index_hd = grp.create_dataset("index", index.shape, index.dtype)

        index_hd.write_direct(index)
        
        columns = np.array(list(data.columns))
        
        columns_hd = grp.create_dataset("columns", columns.shape, columns.dtype)
        
        columns_hd.write_direct(columns)

class GamExperiment(object):
    """A class for storing and processing data associated with a GAM Experiment"""
    
    def __init__(self, hdf5_path, num_processes=1):
        """Open a saved gam_experiment"""
        
        # Open the hdf5 store
        self.store = self.open_store(hdf5_path)
        
        # Get the experimental data from the store
        self.experimental_data = GamExperimentalData(self.store)
        
        # Get the frequency matrix
        self.freq_matrix = GamFrequencyMatrix(self.store)

        # Store the number of processes available for matrix computation
        self.num_processes = num_processes
    
    def get_chrom_processed_matrix(self, chrom, method=None):
        
        starting_index = ( chrom, self.experimental_data.data[chrom].columns[0] )
        stopping_index = ( chrom, self.experimental_data.data[chrom].columns[-1] )
        
        return self.get_index_processed_matrix(starting_index, stopping_index, starting_index, stopping_index, method)
        
    def get_index_processed_matrix(self, index1_start, index1_stop, index2_start, index2_stop, method=None):
        
        loc1_start = self.experimental_data.data.columns.get_loc(index1_start)
        loc1_stop = self.experimental_data.data.columns.get_loc(index1_stop)
        loc2_start = self.experimental_data.data.columns.get_loc(index2_start)
        loc2_stop = self.experimental_data.data.columns.get_loc(index2_stop)
    
        return self.get_loc_processed_matrix(loc1_start, loc1_stop, loc2_start, loc2_stop, method)
    
    def calculate_loc_frequency_matrix(self, loc1_start, loc1_stop, loc2_start, loc2_stop):
        
        index1_list = self.experimental_data.data.columns[loc1_start:loc1_stop]
        index2_list = self.experimental_data.data.columns[loc2_start:loc2_stop]

        data_1 = np.array(self.experimental_data.data.iloc[:,loc1_start:loc1_stop])
        len_1 = len(data_1[0])
        data_2 = np.array(self.experimental_data.data.iloc[:,loc2_start:loc2_stop])
        len_2 = len(data_2[0])
        full_data = np.concatenate((data_1, data_2), axis=1)
        
        combinations = itertools.product(range(len_1), range(len_1, len_1 + len_2))

        from multiprocessing import Manager, Queue, Process
        import time

        print 'Using {0} processes'.format(self.num_processes)
        q = Queue(maxsize=self.num_processes)
        manager = Manager()

        d = manager.dict()
        process_pool = []
        for n in range(self.num_processes):
            process_pool.append(Process(target=f, args=(q,d)))
        
        map(lambda p: p.start(), process_pool)
        
        start_time = time.time()
        for key in combinations:
            loc1, loc2 = key
            data = full_data[:,[loc1, loc2]]
            q.put((key, data))
            
        for i in range(self.num_processes):
            q.put((None, None))
            
        map(lambda p: p.join(), process_pool)

        print 'Time elapsed: ', time.time() - start_time

        freqs = np.zeros((len_1, len_2, 2, 2))

        for key in d.keys():
            loc1, loc2 = key
            loc2 = loc2 - len_1
            freqs[loc1][loc2] = d[key]

        return freqs
        
    def get_loc_frequency_matrix(self, loc1_start, loc1_stop, loc2_start, loc2_stop):
        
        freqs = self.freq_matrix.get_matrix(loc1_start, loc1_stop, loc2_start, loc2_stop)
        
        if freqs is None:
            
            freqs = self.calculate_loc_frequency_matrix(loc1_start, loc1_stop, loc2_start, loc2_stop)
            
            self.freq_matrix.cache_freqs(loc1_start, loc1_stop, loc2_start, loc2_stop, freqs)
            
        return freqs
        
    def get_loc_processed_matrix(self, loc1_start, loc1_stop, loc2_start, loc2_stop, method=None):
        
        if method is None:
            method = self.odds_ratio
        
        index1_list = self.experimental_data.data.columns[loc1_start:loc1_stop]
        index2_list = self.experimental_data.data.columns[loc2_start:loc2_stop]
        
        freqs = self.get_loc_frequency_matrix(loc1_start, loc1_stop, loc2_start, loc2_stop)
        
        stored_shape = freqs.shape[:2]
        
        processed = map(method, freqs.reshape((stored_shape[0] * stored_shape[1],2,2)))
        
        return np.array(processed).reshape(stored_shape)
    
    def odds_ratio(self, N):
        return np.log(N[0,0]) + np.log(N[1,1]) - np.log(np.log(N[0,1])) - np.log(N[1,0])
        
    @staticmethod
    def from_multibam(segmentation_multibam, hdf5_path, num_processes=1):
        """Create a new experiment from a segmentation multibam file"""
        
        # Create the hdf5 store
        store = GamExperiment.create_store(hdf5_path)
        
        # Create the experimental data in the datastore
        experimental_data = GamExperimentalData.from_multibam(segmentation_multibam, store)
        
        # Create a new frequency matrix
        freq_matrix = GamFrequencyMatrix.from_no_windows(store, experimental_data.no_windows)
        
        # close the store
        store.close()
        
        # Return the object
        return GamExperiment(hdf5_path, num_processes)
    
    @staticmethod
    def create_store(hdf5_path):
        """Create a new hdf5 store without overwriting an old one"""
        
        if os.path.exists(hdf5_path):
            # Don't overwrite an existing hdf5 file
            raise HDF5FileExistsError('The file {} already exists and would be overwritten by this operation. Please provide the path for a new HDF5 file'.format(hdf5_path))
        
        # Create the store
        store = h5py.File(hdf5_path,'w')
        
        return store
    
    def open_store(self, hdf5_path):
        """Open an hdf5 store and return the store object"""
        
        if not os.path.exists(hdf5_path):
            raise IOError("[Errno 2] No such file or directory: '{}'".format(hdf5_path))
            
        store = h5py.File(hdf5_path,'r+')
        
        return store
    
    def close(self):
        
        self.store.close()
        
def count_frequency(samples):
    """Take a table of two columnds and return [[ no_both_present, no_1_only],[ no_2_only, no_neither_present]]"""
    
    counts = np.array([[0,0], [0,0]])
    
    for s in samples:

        counts[s[0]][s[1]] += 1
        
    return counts 

def f(input_queue,result_dict):
    key, data = input_queue.get()
    while key is not None:
        result_dict[key] = count_frequency(data)
        key,data = input_queue.get()

