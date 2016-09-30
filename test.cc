#include <Sequence/PolyTable.hpp>
#include <iostream>
#include <algorithm>
#include <iterator>

using namespace Sequence;
using namespace std;


int
main(int argc, char **argv)
{
    vector<pair<double, string>> data;

    data.emplace_back(0.1, "01001");
    data.emplace_back(0.2, "10110");

    PolyTable p(data);

    for (unsigned i = 0; i < p.size(); ++i)
        {
            auto x = p[i];
            //auto v = &x.vector;

            for (std::size_t j = 0; j < x.size(); ++j)
                {
                    std::cout << x[j] << ' ';
                }
            std::cout << '\n';
			for_each(x.begin(),x.end(),[](const char c){ std::cout << c << ' '; });
			std::cout << '\n';
        }

    for (std::size_t i = 0; i < p.numsites(); ++i)
        {
            auto x = p.site(i);
			std::cout << typeid(x).name() << '\n';
            //auto v = &x.vector;

            for (std::size_t j = 0; j < x.size(); ++j)
                {
                    std::cout << x[j]  << ' ';
                }
            std::cout << '\n';
			auto y = x.begin();
			std::cout << typeid(decltype(y)::value_type()).name() << '\n';
			std::cout << std::iterator_traits<decltype(y)>::value_type() << '\n';
			for_each(x.begin(),x.end(),[](const char c){ std::cout << c << ' '; });
			std::cout << '\n';
			std::copy(x.begin(),x.end(),ostream_iterator<char>(std::cout," -> "));
			std::cout << '\n';
			std::string s(x.begin(),x.end());
			std::cout << s << '\n';
        }
}
