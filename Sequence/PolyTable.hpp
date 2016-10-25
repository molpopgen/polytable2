#ifndef SEQUENCE_POLYTABLE_HPP_
#define SEQUENCE_POLYTABLE_HPP_

#include <gsl/gsl_matrix_char.h>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
namespace Sequence
{
    struct delete_matrix
    {
        inline void
        operator()(gsl_matrix_char *m) const noexcept
        {
            gsl_matrix_char_free(m);
        }
    };

    template <typename T> struct gsl_vector_view_wrapper
    /*!
     * \brief Iterable wrapper around GSL vector views.
     */
    {
        using element_type = T;
        T view;
        gsl_vector_view_wrapper(T t) : view(std::move(t)) {}

        using value_type =
            typename std::remove_pointer<decltype(view.vector.data)>::type;

        struct iterator_wrapper
        {
            decltype(view.vector.data) x;
            std::size_t stride;
            using value_type =
                typename std::iterator_traits<decltype(x)>::value_type;
            using reference = value_type &;
            using pointer = value_type *;
            using difference_type =
                typename std::iterator_traits<pointer>::difference_type;
            using iterator_category = std::random_access_iterator_tag;
            inline value_type operator*() { return *x; }
            inline value_type operator*() const { return *x; }
            iterator_wrapper(pointer t, std::size_t stride_)
                : x(t), stride(stride_)
            {
            }
            iterator_wrapper &operator++()
            {
                x += stride;
                return *this;
            }
            iterator_wrapper &operator--()
            {
                x -= stride;
                return *this;
            }
            iterator_wrapper &
            operator+(difference_type d)
            {
                return iterator_wrapper(*this) += d;
            }
            iterator_wrapper &
            operator+=(difference_type d)
            {
                this->x += d;
                return *this;
            }
            iterator_wrapper &
            operator-(difference_type d)
            {
                return iterator_wrapper(*this) -= d;
            }
            iterator_wrapper &
            operator-=(difference_type d)
            {
                this->x -= d;
                return *this;
            }
            difference_type
            operator-(iterator_wrapper itr) const
            {
                return this->x - itr.x;
            }
            pointer operator->() { return x; }
            pointer operator->() const { return x; }

            bool
            operator!=(const iterator_wrapper &rhs)
            {
                return this->x != rhs.x;
            }
            bool
            operator==(const iterator_wrapper &rhs)
            {
                return this->x == rhs.x;
            }
            bool
            operator<(const iterator_wrapper &rhs)
            {
                return this->x < rhs.x;
            }
            bool
            operator>(const iterator_wrapper &rhs)
            {
                return this->x > rhs.x;
            }
        };
        using iterator = iterator_wrapper;

        using const_iterator = typename std::add_const<iterator>::type;

        std::size_t
        size() const
        {
            return view.vector.size;
        }

        value_type operator[](const std::size_t i) const
        {
            return gsl_vector_char_get(&view.vector, i);
        }
        iterator
        begin()
        {
            return iterator(view.vector.data, view.vector.stride);
        }
        iterator
        end()
        {
            return iterator(view.vector.data
                                + view.vector.size * view.vector.stride,
                            view.vector.stride);
        }
        const_iterator
        begin() const
        {
            return const_iterator(view.vector.data, view.vector.stride);
        }
        const_iterator
        end() const
        {
            return const_iterator(view.vector.data
                                      + view.vector.size * view.vector.stride,
                                  view.vector.stride);
        }
        const_iterator
        cbegin() const
        {
            return view.vector.data;
        }
        const_iterator
        cend() const
        {
            return view.vector.data + view.vector.size * view.vector.stride;
        }
    };

    class PolyTable
    {
      private:
        using pimpl_t = std::unique_ptr<gsl_matrix_char, delete_matrix>;
        pimpl_t impl;
        std::vector<double> positions;

      public:
        PolyTable(const std::vector<std::pair<double, std::string>> &sites);
        //! Alias for data view wrapper
        using view_wrapper = gsl_vector_view_wrapper<gsl_vector_char_view>;
        static_assert(
            !std::is_const<typename view_wrapper::element_type>::value,
            "view_wrapper::element_type must be const");
        //! Alias for const data view wrapper
        using const_view_wrapper
            = gsl_vector_view_wrapper<gsl_vector_char_const_view>;
        static_assert(
            std::is_const<typename const_view_wrapper::element_type>::value,
            "const_view_wrapper::element_type must be const");
        //! Non-const view for a haplotype/sequence
        using haplotype_view = view_wrapper;
        //! Const view for a haplotype/sequence
        using const_haplotype_view = const_view_wrapper;
        //! Non-const view of a column/variable site
        using site_view = std::pair<double, view_wrapper>;
        //! Const view of a column/variable site
        using const_site_view = std::pair<const double, const_view_wrapper>;
        //! \return Number of columns (mutation positions)
        std::size_t numsites() const;
        //! \return Number of rows (sequences)
        std::size_t size() const;

        // Member access

        /*!
         * \returns the i-th sequence as a PolyTable::haplotype_view
         * \warning Not range-checked
         */
        haplotype_view operator[](const std::size_t i);
        /*!
         * \returns the i-th sequence as a PolyTable::const_haplotype_view
         * \warning Not range-checked
         */
        const_haplotype_view operator[](const std::size_t i) const;

		//! Range-checked access to i-th PolyTable::haplotype_view
		haplotype_view at(const std::size_t i);
		//! Range-checked access to i-th PolyTable::const_haplotype_view
		const_haplotype_view at(const std::size_t i) const;
        /*!
         * \returns the i-th segregating site as a PolyTable::site_view
         * \warning Not range-checked
         */
        site_view site(const std::size_t i);
        /*!
         * \returns the i-th segregating site as a PolyTable::const_site_view
         * \warning Not range-checked
         */
        const_site_view site(const std::size_t i) const;
    };
}

#endif
