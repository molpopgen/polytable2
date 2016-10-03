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
	{
		using element_type = T;
        T view;
        gsl_vector_view_wrapper(T t) : view(std::move(t)) {}

		using value_type = typename std::remove_pointer<decltype(view.vector.data)>::type;

		struct iterator_wrapper
        {
            decltype(view.vector.data) x;
            std::size_t stride;
            using value_type = typename std::iterator_traits<decltype(x)>::value_type;
			using reference = value_type &;
			using pointer = value_type *;
			using difference_type = typename std::iterator_traits<pointer>::difference_type;
			//using iterator_category = std::random_access_iterator_tag;
			using iterator_category = std::forward_iterator_tag;
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
            bool
            operator!=(const iterator_wrapper &rhs)
            {
                return this->x != rhs.x;
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
		using view_wrapper = gsl_vector_view_wrapper<gsl_vector_char_view>;
		static_assert( !std::is_const<typename view_wrapper::element_type>::value,
				"view_wrapper::element_type must be const");
		using const_view_wrapper 
            = gsl_vector_view_wrapper<gsl_vector_char_const_view>;
		static_assert( std::is_const<typename const_view_wrapper::element_type>::value,
				"const_view_wrapper::element_type must be const");
        using haplotype_view = view_wrapper;
        using const_haplotype_view = const_view_wrapper;
        using site_view = std::pair<double,view_wrapper>;
        using const_site_view
            = std::pair<double,const_view_wrapper>;
        std::size_t numsites() const;
        std::size_t size() const;

        // Member access

        //! Return the i-th haplotype
        haplotype_view operator[](const std::size_t i);
        const_haplotype_view operator[](const std::size_t i) const;

        site_view site(const std::size_t i);
        const_site_view site(const std::size_t i) const;
    };
}

#endif
