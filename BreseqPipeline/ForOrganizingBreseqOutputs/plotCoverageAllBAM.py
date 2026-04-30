import os
import pysam
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.io # pip install scipy

def process_bam_file(strFileName):
    # Start timer
    start_time = time.time()

    # Load BAM file
    bamFile = pysam.AlignmentFile(strFileName + '.bam', 'rb')
    firstRead = next(bamFile)

    print(f"Header: {bamFile.header}")
    print(f"References: {bamFile.references}")
    print(f"Lengths: {bamFile.lengths}")

    # Define reference sequence name
    refName = bamFile.references[0]

    # Get reference sequence length
    refLength = bamFile.lengths[bamFile.references.index(refName)]

    # Define moving average window
    windowSize = 10000

    # Initialize arrays for coverage and position
    coverage = np.zeros(refLength, dtype=int)
    position = np.arange(refLength)
    read_counter = 0 

    firstRead = next(bamFile)

    # Loop through each read in BAM file
    for read in bamFile:
        # Skip unmapped reads
        if read.is_unmapped:
            continue

        # Skip reads with mapping quality less than 10
        if read.mapping_quality < 10:
            continue

        # Get read start and end positions
        start = read.reference_start
        end = read.reference_end

        # Ensure that read alignment positions are within the bounds of the reference sequence
        start = max(start, 0)
        end = min(end, refLength)

        # Add coverage to each position in the window
        for j in range(start, end):
            coverage[j] += 1

        # Increment read counter
        read_counter += 1

        if read_counter % 1000000 == 0:
            print(f"Processed {read_counter/1000000}M reads")

        # Stop after reading 10000 reads (for debug)
        #if read_counter >= 10000:
        #   break

    # Smooth coverage using moving average
    smoothCoverage = np.convolve(coverage, np.ones(windowSize) / windowSize, mode='same')

    # Save coverage as MATLAB-readable file
    scipy.io.savemat(strFileName + '_coverage.mat', {'coverage': coverage})

    # Close BAM file
    bamFile.close()

    # End timer and print elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")

    ## Plot coverage
    #plt.plot(position, smoothCoverage)
    #plt.xlabel('Position')
    #plt.ylabel('Coverage')
    #plt.title('Coverage Across DNA Sequence')
    #plt.show()

def main():
    # Get all BAM files in the current directory
    bam_files = [f for f in os.listdir('.') if f.endswith('.bam')]

    # Process each BAM file
    for bam_file in bam_files:
        # Get file name without extension
        strFileName = os.path.splitext(bam_file)[0]

        # Process BAM file
        process_bam_file(strFileName)

if __name__ == '__main__':
    main()
