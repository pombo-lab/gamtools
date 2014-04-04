import numpy as np
import os
import itertools
import h5py
from multiprocessing import Pool
from bisect import bisect_left


def D(n):
    total = n.sum()
    f = n / float(total)
    p1 = f[0][1] + f[1][1]
    p2 = f[0][0] + f[1][0]
    q1 = f[1][0] + f[1][1]
    q2 = f[0][0] + f[0][1]
    expected = p1 * q1
    return f[1][1] - expected


def Dmax(n):
    total = n.sum()
    f = n / float(total)
    p1 = f[0][1] + f[1][1]
    p2 = f[0][0] + f[1][0]
    q1 = f[1][0] + f[1][1]
    q2 = f[0][0] + f[0][1]
    d = D(n)
    if d > 0:
        return min([p1 * q2, p2 * q1])
    elif d < 0:
        return min([p1 * q1, p2 * q2])


def Dprime(n):
    d = D(n)
    dmax = Dmax(n)
    if dmax is None:
        return 0.0
    else:
        return d / dmax


def corr(n):
    d = D(n)
    total = n.sum()
    f = n / float(total)
    p1 = f[0][1] + f[1][1]
    p2 = f[0][0] + f[1][0]
    q1 = f[1][0] + f[1][1]
    q2 = f[0][0] + f[0][1]
    return d / np.sqrt(p1 * p2 * q1 * q2)


class HDF5FileExistsError(Exception):
    pass


class Hdf5StoreRequiresCacheException(Exception):
    pass


class GamFrequencyMatrix(object):

    """A class for abstracting access to the HDF5 store underlying a frequency matrix"""

    def __init__(self, data):

        self.data = data

    @staticmethod
    def from_store(hdf5_store):
        """Get a frequency matrix from the store"""

        data = hdf5_store.store["processed_data"]["frequencies"]
        return GamFrequencyMatrix(data)

    @staticmethod
    def from_no_windows(hdf5_store, no_windows):

        data = hdf5_store.create_freq_matrix(no_windows)

        return GamFrequencyMatrix(data)

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


class MultibamProcessor(object):

    """A class for extracting GAM experimental data from the multibam style segmentation file."""

    def __init__(self):
        super(MultibamProcessor, self).__init__()

    def read(self, multibam_pointer):
        """Actually read data from the multibam pointer"""

        windows = []
        data = []
        for lnum, line in enumerate(multibam_pointer):
            if lnum == 0:
                sample_names = self.extract_sample_names(line)
            else:
                window, segmentation = self.parse_line(line)
                windows.append(window)
                data.append(segmentation)

        data = np.array(data).transpose()
        data = np.ascontiguousarray(data, data.dtype)

        return data, sample_names, windows

    def parse_line(self, line):
        """Extract information from one line of the multibam"""

        fields = line.strip().split()
        chrom, start, end = fields[:3]
        start, end = int(start), int(end)
        window = (chrom, start, end)
        segmentation = map(int, fields[3:])

        return window, segmentation

    def extract_sample_name(self, file_path):
        basename = os.path.basename(file_path)
        return basename.split('.')[0]

    def extract_sample_names(self, header_line):
        return map(self.extract_sample_name, header_line.split()[2:])


class GamHdf5Store(object):

    def __init__(self, store, compression):

        self.store, self.compression = store, compression

    @staticmethod
    def create_store(hdf5_path, compression=None):
        """Create a new hdf5 store without overwriting an old one"""

        if os.path.exists(hdf5_path):
            # Don't overwrite an existing hdf5 file
            raise HDF5FileExistsError(
                'The file {0} already exists and would be overwritten by this operation. Please provide the path for a new HDF5 file'.format(hdf5_path))

        # Create the store
        store = h5py.File(hdf5_path, 'w')

        return GamHdf5Store(store, compression)

    @staticmethod
    def open_store(hdf5_path):
        """Open an hdf5 store and return the store object"""

        if not os.path.exists(hdf5_path):
            raise IOError(
                "[Errno 2] No such file or directory: '{0}'".format(hdf5_path))

        store = h5py.File(hdf5_path, 'r+')

        # Get the compression
        compression = store['processed_data']['frequencies'].compression

        return GamHdf5Store(store, compression)

    def data_to_store(self, data, columns, windows):
        """Save experimental data to an hdf5 store"""

        grp = self.store.create_group("experimental_data")

        segmentation_hd = grp.create_dataset(
            "segmentation", data.shape, data.dtype, compression=self.compression)

        segmentation_hd.write_direct(data)

        windows = np.array(windows)

        index_hd = grp.create_dataset("windows", windows.shape, windows.dtype)

        index_hd.write_direct(windows)

        columns = np.array(columns)

        columns_hd = grp.create_dataset(
            "columns", columns.shape, columns.dtype)

        columns_hd.write_direct(columns)

    def convert_window(self, window):

        chrom, start, stop = window

        return chrom, int(start), int(stop)

    def data_from_store(self):
        """Retrieve data from an hdf5 store and use it to recreate the DataFrame"""

        data = np.array(self.store["experimental_data"]["segmentation"][:])
        columns = np.array(self.store["experimental_data"]["columns"][:])
        windows = map(self.convert_window,
                      np.array(self.store["experimental_data"]["windows"][:]))

        return data, columns, windows

    def create_freq_matrix(self, no_windows):
        """Create a frequency matrix dataset stored in the HDF5 file"""

        grp = self.store.create_group("processed_data")

        freq_matrix = grp.create_dataset("frequencies",
                                         (no_windows, no_windows, 2, 2),
                                         dtype='i',
                                         compression=self.compression)

        return freq_matrix


class GamMockHdf5Store(object):

    def data_to_store(self, data, columns, windows):
        """Dont Save experimental data to an hdf5 store"""

        pass

    def create_freq_matrix(self, no_windows):

        class DummyDataStore():

            def __getitem__(*args):

                return np.array([[0]])

            def __setitem__(*args):

                pass

        return DummyDataStore()

    def close(self):

        pass


class GamExperimentalData(object):

    """A class for abstracting access to the original experimental segmentation"""

    def __init__(self, data, columns, windows, hdf5_store):

        self.data, self.columns, self.windows, self.store = data, columns, windows, hdf5_store
        self.no_windows = len(self.windows)

    @staticmethod
    def from_multibam(multibam_path, hdf5_store):

        with open(multibam_path, 'r') as multibam_pointer:
            data, columns, windows = MultibamProcessor().read(multibam_pointer)

        hdf5_store.data_to_store(data, columns, windows)

        return GamExperimentalData(data, columns, windows, hdf5_store)

    @staticmethod
    def from_store(hdf5_store):

        data, columns, windows = hdf5_store.data_from_store()

        return GamExperimentalData(data, columns, windows, hdf5_store)


class GamExperiment(object):

    """A class for storing and processing data associated with a GAM Experiment"""

    def __init__(self, hdf5_store, experimental_data, frequency_data, num_processes=1):
        """Open a saved gam_experiment"""

        # Store the number of processes available for matrix computation
        self.num_processes = num_processes

        # Starting the worker pool forks the process, so do this as early as possible to
        # reduce memory usage
        if self.num_processes > 1:
            self.processor = MultithreadedMatrixProcesser(num_processes)
        else:
            self.processor = MatrixProcesser()

        # Open the hdf5 store
        self.store = hdf5_store

        # Get the experimental data from the store
        self.experimental_data = experimental_data

        # Get the frequency matrix
        self.freq_matrix = frequency_data

    @staticmethod
    def from_multibam(segmentation_multibam, hdf5_path, num_processes=1, compression=None):
        """Create a new experiment from a segmentation multibam file"""

        # Create the hdf5 store
        store = GamHdf5Store.create_store(hdf5_path, compression)

        # Create the experimental data in the datastore
        experimental_data = GamExperimentalData.from_multibam(
            segmentation_multibam, store)

        # Create a new frequency matrix
        freq_matrix = GamFrequencyMatrix.from_no_windows(
            store, experimental_data.no_windows)

        # Return the object
        return GamExperiment(store, experimental_data, freq_matrix, num_processes)

    @staticmethod
    def from_multibam_no_cache(segmentation_multibam, num_processes=1):
        """Create a new experiment from a segmentation multibam file"""

        # Create the hdf5 store
        store = GamMockHdf5Store()

        # Create the experimental data in the datastore
        experimental_data = GamExperimentalData.from_multibam(
            segmentation_multibam, store)

        # Create a new frequency matrix
        freq_matrix = GamFrequencyMatrix.from_no_windows(
            store, experimental_data.no_windows)

        # Return the object
        return GamExperiment(store, experimental_data, freq_matrix, num_processes)

    @staticmethod
    def load(hdf5_path, num_processes=1):

        # Load the hdf5 store
        store = GamHdf5Store.open_store(hdf5_path)

        experimental_data = GamExperimentalData.from_store(store)

        freq_matrix = GamFrequencyMatrix.from_store(store)

        # Return the object
        return GamExperiment(store, experimental_data, freq_matrix, num_processes)

    def calculate_loc_frequency_matrix(self, loc1_start, loc1_stop, loc2_start, loc2_stop):

        data_1 = np.array(self.experimental_data.data[:, loc1_start:loc1_stop])
        len_1 = len(data_1[0])
        data_2 = np.array(self.experimental_data.data[:, loc2_start:loc2_stop])
        len_2 = len(data_2[0])
        full_data = np.concatenate((data_1, data_2), axis=1)

        return self.processor.process(len_1, len_2, full_data)

    def get_loc_frequency_matrix(self, loc1_start, loc1_stop, loc2_start, loc2_stop):

        freqs = self.freq_matrix.get_matrix(
            loc1_start, loc1_stop, loc2_start, loc2_stop)

        if freqs is None:

            freqs = self.calculate_loc_frequency_matrix(
                loc1_start, loc1_stop, loc2_start, loc2_stop)

            self.freq_matrix.cache_freqs(
                loc1_start, loc1_stop, loc2_start, loc2_stop, freqs)

        return freqs

    def get_loc_processed_matrix(self, loc1_start, loc1_stop, loc2_start, loc2_stop, method=None):

        if method is None:
            method = corr

        freqs = self.get_loc_frequency_matrix(
            loc1_start, loc1_stop, loc2_start, loc2_stop)

        stored_shape = freqs.shape[:2]

        processed = map(
            method, freqs.reshape((stored_shape[0] * stored_shape[1], 2, 2)))

        return np.array(processed).reshape(stored_shape)

    def get_chrom_start_stop_bins(self, chrom):

        chrom_bins = map(lambda t: t[0], self.experimental_data.windows)

        start_bin = chrom_bins.index(chrom)
        stop_from_end = chrom_bins[::-1].index(chrom)
        stop_bin = len(chrom_bins) - 1 - stop_from_end

        return start_bin, stop_bin

    def get_bins_from_positions(self, chrom, start_pos, stop_pos):

        chrom_start, chrom_stop = self.get_chrom_start_stop_bins(chrom)
        chrom_windows = self.experimental_data.windows[
            chrom_start:chrom_stop + 1]
        chrom_positions = map(lambda t: t[2], chrom_windows)
        start_bin = bisect_left(chrom_positions, start_pos) + chrom_start
        stop_bin = bisect_left(chrom_positions, stop_pos) + chrom_start

        return start_bin, stop_bin

    def parse_location_string(self, string):

        chrom_fields = string.split(':')

        chrom = chrom_fields[0]

        if len(chrom_fields) == 1:
            start, stop = self.get_chrom_start_stop_bins(chrom)

        else:

            pos_fields = chrom_fields[1].split('-')

            def convert_position(pos):
                clean_pos = pos.replace(',', '')
                return int(clean_pos)

            start_pos, stop_pos = map(convert_position, pos_fields)
            start, stop = self.get_bins_from_positions(
                chrom, start_pos, stop_pos)

        # Always return at least one bin width
        if start == stop:
            stop = start + 1

        return start, stop

    def parse_locations(self, location1, location2=None):

        loc1_start, loc1_stop = self.parse_location_string(location1)

        if location2 is None:
            loc2_start, loc2_stop = loc1_start, loc1_stop

        else:
            loc2_start, loc2_stop = self.parse_location_string(location2)

        return loc1_start, loc1_stop, loc2_start, loc2_stop

    def frequencies(self, location1, location2=None):

        return self.get_loc_frequency_matrix(*self.parse_locations(location1, location2))

    def distances(self, location1, location2=None, method=None):

        return self.get_loc_processed_matrix(*self.parse_locations(location1, location2), method=method)

    def close(self):

        self.store.close()

    def __enter__(self):

        return self

    def __exit__(self, type, value, traceback):

        self.close()


def count_frequency(samples):
    """Take a table of two columnds and return [[ no_both_present, no_1_only],[ no_2_only, no_neither_present]]"""

    counts = np.array([[0, 0], [0, 0]])

    for s in samples:

        counts[s[0]][s[1]] += 1

    return counts


class MatrixProcesser(object):

    def process(self, len_1, len_2, full_data):

        combinations = itertools.product(
            range(len_1), range(len_1, len_1 + len_2))

        self._data = full_data
        result = map(self.get_frequency, combinations)

        freqs = np.array(result).reshape((len_1, len_2, 2, 2))

        return freqs

    def get_frequency(self, i):

        p, q = i
        return self.count_frequency(self._data[:, [p, q]])

    def count_frequency(self, samples):

        counts = np.array([[0, 0], [0, 0]])

        for s in samples:

            counts[s[0]][s[1]] += 1

        return counts


class MultithreadedMatrixProcesser(MatrixProcesser):

    def __init__(self, num_processes):

        self.num_processes = num_processes

    def process(self, len_1, len_2, full_data):

        combinations = itertools.product(
            range(len_1), range(len_1, len_1 + len_2))

        print 'Using {0} processes'.format(self.num_processes)

        self._data = full_data
        p = Pool(self.num_processes)
        result = p.map(self.get_frequency, combinations)
        p.close()
        p.join()

        print 'Done processing'

        freqs = np.array(result).reshape((len_1, len_2, 2, 2))

        return freqs

# Magic code that I don't understand, copied from
# http://bytes.com/topic/python/answers/552476-why-cant-you-pickle-instancemethods


def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)


def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

import copy_reg
import types
copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

# End magic code
