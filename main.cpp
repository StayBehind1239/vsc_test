#include <algorithm>
#include <iomanip>
#ifndef __GNUC__
#include <ios>
#endif
#include <iostream>
#include <string>
#include <vector>
#include "dice.cpp"

using std::cin;             using std::sort;
using std::cout;            using std::streamsize;
using std::endl;            using std::string;
using std::setprecision;    using std::vector;

int main() {

  int i = 0;

  while(i < 1000000)
  {
    cout << "We are at: " << i << "\n";

    i = roll(i);

  }
  return 0;
}
