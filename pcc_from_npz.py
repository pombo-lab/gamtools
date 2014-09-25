import itertools
import numpy as np
import argparse
import glob
import logging

parser = argparse.ArgumentParser(description='Calculate PCC between different chromosomes stored as npz files.')
parser.add_argument('-p1','--paths1', metavar='FIRST_PATH', required=True, nargs='+', help='One or more input npz files.')
parser.add_argument('-p2','--paths2', metavar='SECOND_PATH', required=True, nargs='+', help='One or more input npz files.')
parser.add_argument('-d','--distance', metavar='DISTANCE', type=int, default=np.Inf, help='Only consider bins within this distance of each other.')
parser.add_argument('--debug',
    help='Print lots of debugging statements',
    action="store_const",dest="loglevel",const=logging.DEBUG,
    default=logging.WARNING
)
parser.add_argument('--verbose',
    help='Be verbose',
    action="store_const",dest="loglevel",const=logging.INFO
)


class ZipChroms(object):
    
    def __init__(self, path1, path2, distance=np.Inf):
        
        logging.debug( 'initializing')
        self.distance = distance
        logging.debug( 'distance to consider is {0}'.format(self.distance))
        self.path1, self.path2 = path1, path2
        self.loaded = False
        self.array1, self.array2 = None, None
        self.p, self.q = 0, 0
        
    def __iter__(self):
        if not self.loaded:
            logging.debug( 'loading data')
            self.load_matrices()
        return self
        
    def next(self):
           
        while True:
            try:
                result = self.get_next()
            except StopIteration:
                logging.debug( 'Deleting')
                del self.array1
                del self.array2
                self.loaded = False
                raise StopIteration
            if result:
                return result
            
    def get_next(self):
        
        try:
            xi = self.array1.next()
            yi = self.array2.next()
        except StopIteration:
            raise StopIteration
        
        if np.isfinite(xi) and np.isfinite(yi):
            return xi, yi
        else:
            return None
        
    def load_matrices(self):
        
        self.array1 = self.load_matrix(self.path1)
        self.array2 = self.load_matrix(self.path2)
        self.loaded = True
        
    def load_matrix(self, path):
        
        arr = np.load(path)['scores']
        old_shape = arr.shape
        indices = itertools.product(range(old_shape[0]), range(old_shape[1]))
        
        return ( arr[p,q] for p,q in indices if np.abs((p - q) < self.distance) )

def get_means(iterable):
    
    n, x_sum, y_sum = 0, 0, 0
    for x,y in iterable:
        n += 1
        x_sum += x
        y_sum += y
    return n, x_sum / n, y_sum / n

def get_sdevs(iterable, n, x_mean, y_mean):
    
    x_ssqrs, y_ssqrs = 0, 0
    
    for x, y in iterable:
        x_ssqrs += (x - x_mean)**2
        y_ssqrs += (y - y_mean)**2
        
    x_var = x_ssqrs / (n-1)
    y_var = y_ssqrs / (n-1)
    
    return np.sqrt(x_var), np.sqrt(y_var)

def do_pcc(iterable, n, x_mean, y_mean, x_sd, y_sd):
    
    prod = 0
    
    for x, y in iterable:
        x_ss = (x - x_mean) / x_sd
        y_ss = (y - y_mean) / y_sd
        prod += x_ss * y_ss
        
    return prod / (n-1)

def get_pcc(paths1, paths2, distance=np.Inf):
    
    _zipped_paths = zip(paths1, paths2)
    
    _zipped_arrays = [ ZipChroms(p1, p2, distance=distance) for p1, p2 in _zipped_paths ]
    
    n, x_mean, y_mean = get_means(itertools.chain(*_zipped_arrays))
    
    x_sd, y_sd = get_sdevs(itertools.chain(*_zipped_arrays), n, x_mean, y_mean)
    
    r = do_pcc(itertools.chain(*_zipped_arrays), n, x_mean, y_mean, x_sd, y_sd)
    
    return r


def main(args):

    if len(args.paths1) > 1:
        args.paths1.sort()
        args.paths2.sort()
    else:
        args.paths1 = sorted(glob.glob(args.paths1[0]))
        args.paths2 = sorted(glob.glob(args.paths2[0]))

    for paths in zip(args.paths1, args.paths2):
        logging.info('{0} <--> {1}'.format(*paths))

    r = get_pcc(args.paths1, args.paths2, distance=args.distance)

    print 'PCC is: {0}'.format(r)

if __name__ == '__main__':

    args = parser.parse_args()

    logging.basicConfig(level=args.loglevel)

    main(args)

