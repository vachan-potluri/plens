/**
 * @file metric_terms.h
 * @brief A class for storing DGSEM metric terms of a cell
 */

#include "metric_terms.h"

/**
 * Constructor using `FEValues`. Calls reinit() internally
 */
template <int dim>
MetricTerms<dim>::MetricTerms(const FEValues<dim>& fev)
{
    reinit(fev);
}



/**
 * The main function: computes and stores all the metric terms according to `fev` provided.
 */
template <int dim>
void MetricTerms<dim>::reinit(const FEValues<dim>& fev)
{}
