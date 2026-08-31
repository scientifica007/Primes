#include <primesieve.hpp>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>

int main(int argc, char** argv)
{
    std::uint64_t limit = 100000000ULL;
    if (argc >= 2)
        limit = std::strtoull(argv[1], nullptr, 10);

    auto start = std::chrono::steady_clock::now();

    primesieve::iterator it;
    std::uint64_t count = 0;
    std::uint64_t checksum = 0;
    std::uint64_t last = 0;

    for (std::uint64_t p = it.next_prime(); p <= limit; p = it.next_prime()) {
        ++count;
        checksum += p;
        last = p;
    }

    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;

    std::cout << "limit=" << limit
              << " count=" << count
              << " checksum=" << checksum
              << " last=" << last
              << " seconds=" << std::fixed << std::setprecision(9) << elapsed.count()
              << '\n';
    return 0;
}
