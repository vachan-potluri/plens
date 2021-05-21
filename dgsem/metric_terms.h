/**
 * @file metric_terms.h
 * @brief A class for storing DGSEM metric terms of a cell
 */

#ifndef METRIC_TERMS_H
#define METRIC_TERMS_H

/**
 * @class MetricTerms
 * This class computes and stores metric terms required for DGSEM algorithm for a cell. Using this
 * class is simple. First, construct the class by specifying the polynomial degree. Then, call
 * MetricTerms::reinit() by passing an appropriate `FEValues` object. When making this call, the
 * `FEValues` object must also have been re-initialised on the desired cell.
 *
 * The most common way of using this with a `DoFHandler` is to use a `std::map` of `MetricTerms`
 * objects using the cell index as the key. This will however not work for adaptively refined
 * meshes where a cell index is not sufficient to uniquely identify a cell. However, we are not
 * interested in such cases and for our purposes, a simple map will suffice.
 */

#endif
