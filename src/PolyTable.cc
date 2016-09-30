#include <Sequence/PolyTable.hpp>

#include <algorithm>
#include <iostream>
namespace Sequence
{
    PolyTable::PolyTable(
        const std::vector<std::pair<double, std::string>> &sites)
        : impl(nullptr)
    {
        if (!sites.empty())
            {
                impl.reset(gsl_matrix_char_alloc(sites.size(),
                                                 sites[0].second.size()));
            }
        std::size_t offset = 0;
        for (auto &&si : sites)
            {
                std::copy(si.second.begin(), si.second.end(),
                          impl->data + offset);
                offset += si.second.size();
            }
        for (std::size_t i = 0; i < impl->size1; ++i)
            {

                for (std::size_t j = 0; j < impl->size2; ++j)
                    {
                        std::cout << gsl_matrix_char_get(impl.get(), i, j)
                                  << ' ';
                    }
                std::cout << '\n';
            }
    }

    std::size_t
    PolyTable::size() const
    {
        return impl->size2;
    }

    std::size_t
    PolyTable::numsites() const
    {
        return impl->size1;
    }

	PolyTable::view_type PolyTable::operator[](const std::size_t i)
    {
        return PolyTable::view_type(gsl_matrix_char_column(impl.get(), i));
    }

	PolyTable::const_view_type PolyTable::operator[](const std::size_t i) const
    {
        return PolyTable::const_view_type(gsl_matrix_char_const_column(impl.get(), i));
    }

	PolyTable::view_type
    PolyTable::site(const std::size_t i)
    {
        return PolyTable::view_type(gsl_matrix_char_row(impl.get(), i));
    }

	PolyTable::const_view_type
    PolyTable::site(const std::size_t i) const
    {
        return PolyTable::const_view_type(gsl_matrix_char_const_row(impl.get(), i));
    }
}
