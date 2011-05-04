#ifndef PEAK_H
#define PEAK_H

/**
 * AUTHOR: Kha Nguyen
 * CREATE DATE: 04/21/2011
 *
 * Object for representing one peak in a spectrum
 *
 * A peak is primarily identified via its intensity (height) and location \
 * (position on the m/z axis)
 *
 * \file Peak.h
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <vector>
#include "objects.h"

class Peak {
public:
    /**
     * Construct a Peak object
     * @param intensity: intensity for the new peak
     * @param location: location for the new peak
     */
    Peak(FLOAT_T intensity, FLOAT_T location);
    
    /**
     * Return the intensity of this Peak
     */
    FLOAT_T getIntensity();
    
    /**
     * Return the intensity rank of this Peak
     */
    FLOAT_T getIntensityRank();
    
    /**
     * Return the location of this Peak
     */
    FLOAT_T getLocation();
    
    /**
     * Set the intensity of this Peak
     */
    void setIntensity(FLOAT_T intensity);
    
    /**
     * Set the intensity rank of this Peak
     */
    void setIntensityRank(FLOAT_T intensity_rank);
    
    /**
     * Set the location of this Peak
     */
    void setLocation(FLOAT_T location);
    
    /**
     * Print the intensity and location of this peak to stdout
     */
    void print();
    
    /**
     * Compare the intensity of this Peak and another Peak
     * Return true if this Peak is greater, false otherwise
     */
    bool compareByIntensity(Peak other);
    
    /**
     * Compare the mz(location) of this Peak and another Peak
     * Return true if the other Peak is greater, false otherwise
     */
    bool compareByMZ(Peak other);
    
private:
    FLOAT_T intensity_;          // The intensity of this peak
    FLOAT_T intensity_rank_;     // The rank intensity of this peak
    FLOAT_T location_;           // The location of this peak
    
};

/**
 * Return a vector of allocated Peak pointers
 */
std::vector<Peak*> allocate_peak_vector(unsigned int num_peaks);

/**
 * Free a Peak pointer array
 */
void free_peak_vector(std::vector<Peak*> &peaks);

/**
 * Sort peaks by their intensity or location
 * Use the lib function sort()
 */
void sort_peaks(std::vector<Peak*> &peak_array, PEAK_SORT_TYPE_T sort_type);

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

#endif
