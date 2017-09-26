#include "readToVec.hpp"

std::vector< std::vector<double> >  readToVec(std::string filename)
{
    std::vector< std::vector<double> >  values;
    std::ifstream fin(filename.c_str());
    for (std::string line; std::getline(fin, line); )
    {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream in(line);
        values.push_back(
            std::vector<double>(std::istream_iterator<double>(in),
                                std::istream_iterator<double>()));
    }

    for (std::vector< std::vector<double> > ::const_iterator
             it(values.begin()), end(values.end()); it != end; ++it) {
        std::copy(it->begin(), it->end(),
                  std::ostream_iterator<double>(std::cout, ", "));
        std::cout << "\n";
    }
    return values;
}
