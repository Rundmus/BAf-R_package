# BAf-R_package
BAf - an R package to handle the data obtained from binder arrays in Plasma profiling facility

This package has been developed to facilitate other collaborators analysis and minimize errors in data handling. The data from binder arrays are composed of heterogenous informations that cannot be stored in a data frame. They are, for examples, sample information, signal values, binder information, and assay specific data. When we store those data in a few data frames and want to extract a subset of the data, it is prone to make a mistake, e.g. forgetting to get a subset of assay specific data. In order to minimize such errors, this package provides classes and methods, which let users avoid such errors. 
