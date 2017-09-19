#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#ifndef READTOVEC_HPP_INCLUDED
#define READTOVEC_HPP_INCLUDED

/// read a n x m excel-style csv file to a vector of vectors (containing doubles, easily changed to other types)
///     input: file -> string naming the file to be read in
std::vector< std::vector<double> > readToVec(std::string file);

#endif // READTOVEC_HPP_INCLUDED
