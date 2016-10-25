#include <Sequence/PolyTable.hpp>

#include <algorithm>
#include <iostream>
namespace Sequence
{
    PolyTable::PolyTable(
        const std::vector<std::pair<double, std::string>> &sites)
        : impl(nullptr), positions(std::vector<double>())
    {
        if (!sites.empty())
            {
                impl.reset(gsl_matrix_char_alloc(sites.size(),
                                                 sites[0].second.size()));
            }
        std::size_t offset = 0;
        for (auto &&si : sites)
            {
                positions.push_back(si.first);
                std::copy(si.second.begin(), si.second.end(),
                          impl->data + offset);
                offset += si.second.size();
            }
    }

    PolyTable::PolyTable(const PolyTable &t)
        : impl(gsl_matrix_char_alloc(t.size(), t.numsites())),
          positions(t.positions)
    {
        gsl_matrix_char_memcpy(impl.get(), t.impl.get());
    }

    PolyTable::PolyTable(PolyTable &&t)
        : impl(std::move(t.impl)), positions(std::move(t.positions))
    {
        t.impl.reset(nullptr);
    }

    std::size_t
    PolyTable::size() const
    {
        if (impl.get() == nullptr)
            return 0;
        return impl->size2;
    }

    std::size_t
    PolyTable::numsites() const
    {
        if (impl.get() == nullptr)
            return 0;
        return impl->size1;
    }

    PolyTable::haplotype_view PolyTable::operator[](const std::size_t i)
    {
        return PolyTable::haplotype_view(
            gsl_matrix_char_column(impl.get(), i));
    }

    PolyTable::haplotype_view
    PolyTable::at(const std::size_t i)
    {
        if (i >= this->size())
            throw std::out_of_range("PolyTable::at, index out of range");
        return this->operator[](i);
    }

    PolyTable::const_haplotype_view
    PolyTable::at(const std::size_t i) const
    {
        if (i >= this->size())
            throw std::out_of_range("PolyTable::at, index out of range");
        return this->operator[](i);
    }

    PolyTable::const_haplotype_view PolyTable::
    operator[](const std::size_t i) const
    {
        return PolyTable::const_haplotype_view(
            gsl_matrix_char_const_column(impl.get(), i));
    }

    PolyTable::site_view
    PolyTable::site(const std::size_t i)
    {
        return std::make_pair(
            positions[i], view_wrapper(gsl_matrix_char_row(impl.get(), i)));
    }

    PolyTable::const_site_view
    PolyTable::site(const std::size_t i) const
    {
        return std::make_pair(
            positions[i],
            const_view_wrapper(gsl_matrix_char_const_row(impl.get(), i)));
    }

    PolyTable::site_view
    PolyTable::at_site(const std::size_t i)
    {
        if (i >= this->numsites())
            throw std::out_of_range("PolyTable::at_site, index out of range");
        return this->site(i);
    }
    PolyTable::const_site_view
    PolyTable::at_site(const std::size_t i) const
    {
        if (i >= this->numsites())
            throw std::out_of_range("PolyTable::at_site, index out of range");
        return this->site(i);
    }
    void
    PolyTable::swap(PolyTable &t)
    {
        std::swap(this->impl, t.impl);
        this->positions.swap(t.positions);
    }
    void
    swap(PolyTable &a, PolyTable &b)
    {
        a.swap(b);
    }
}

namespace std
{
    template <> void swap(Sequence::PolyTable &a, Sequence::PolyTable &b)
	{
		a.swap(b);
	}
}
