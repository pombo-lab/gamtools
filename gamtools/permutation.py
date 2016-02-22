from . import segmentation
import numpy as np
import pandas as pd

def permute_segmentation(input_segmentation):

    # Only permute positions that were mapped at least once
    # This means we preserve the locations of centromeres etc,
    # but requires that the input data has a large number of
    # columns
    mappable = input_segmentation.sum(axis=1).astype(bool)
    no_windows, no_samples = input_segmentation[mappable].shape

    # Make a copy of the original data
    permutation = input_segmentation.copy()

    # Loop over columns
    for i in range(no_samples):

        # Choose a random position to break the segmentation in two,
        # Swap the two chunks around and write them to the copied df
        offset = np.random.randint(no_windows)

        new_start = input_segmentation[mappable].iloc[offset:,i]
        new_end = input_segmentation[mappable].iloc[:offset,i]
        new_col = list(pd.concat([new_start, new_end]))

        permutation.ix[mappable,i] = new_col

    return permutation

def permute_segmentation_autosomal(input_segmentation):

    # Sex chromosomes, unjoined contigs and mtDNA behave weirdly,
    # so don't permute them into the autosomal regions.

    autosomes = ['chr{0}'.format(c) for c in range(1,20)]

    is_autosomal = input_segmentation.index.get_level_values(0).isin(autosomes)

    segmentation_to_permute = input_segmentation.loc[is_autosomal]

    permuted_segmentation = permute_segmentation(segmentation_to_permute)

    return permuted_segmentation


def permute_segmentation_from_args(args):

    input_segmentation = segmentation.open_segmentation(args.segmentation_file)

    permuted_segmentation = permute_segmentation_autosomal(input_segmentation)

    permuted_segmentation.to_csv(args.output_file, sep='\t', header=True, index=True)
