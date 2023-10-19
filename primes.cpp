#include "CS207/Util.hpp"
#include <math.h>       /* sqrt */
#include <vector>

using namespace std;

/** Return true iff @a n is prime.
 * @pre @a n >= 0 */
bool is_prime(int n)
{
  assert(n >= 0);

  static vector<int> memo = {2}; // vector that stores candidate divisors
  static int max_n_tested = 2;
  
  // increase size of memo if needed
  for (int i = max_n_tested+1; i <= n; ++i)
    memo.push_back(i);
  if ( n > max_n_tested )
    max_n_tested = n;

  vector<int>::const_iterator it;
  int terminating_value = int(sqrt(n)) + 1;
  for (it = memo.begin(); *it < terminating_value; ++it)
    // test division by candidate divisors (which will be prime numbers if is_prime is called for increasing values of n = 2,3,... else it will include additional candidates)
      if (n % *it == 0) {
        // remove n efficiently if it is last element in memo -- it will be if is_prime is called for increasing values of n
        if (memo.back() == n)
          memo.pop_back();
        return false;
      };
  
  if ( n < 2 )
    return false;
  
  return true; 
}


int main()
{
  while (!cin.eof()) {
    // How many primes to test? And should we print them?
    cerr << "Input Number: ";
    int n = 0;
    CS207::getline_parsed(cin, n);
    if (n <= 0)
      break;

    cerr << "Print Primes (y/n): ";
    char confirm = 'n';
    CS207::getline_parsed(cin, confirm);
    bool print_primes = (confirm == 'y' || confirm == 'Y');

    CS207::Clock timer;

    // Loop and count primes from 2 up to n
    int num_primes = 0;
    for (int i = 2; i <= n; ++i) {
      if (is_prime(i)) {
        ++num_primes;
        if (print_primes)
          cout << i << endl;
      }
    }

    double elapsed_time = timer.seconds();

    cout << "There are " << num_primes
         << " primes less than or equal to " << n << ".\n"
         << "Found in " << (1000 * elapsed_time) << " milliseconds.\n\n";
  }

  return 0;
}
