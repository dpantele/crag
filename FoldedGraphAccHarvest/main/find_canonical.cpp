/**
* Simple procedure which finds a canonical pair for the input argumnets
*/

#include <string.h>
#include "acc.h"

using namespace crag;

int main(int argc, char** argv) {
  if (argc != 3) {
    auto name = strrchr(argv[0], '/') + 1;
    if (!name) {
      name = strrchr(argv[0], '\\') + 1;
    }
    if (!name) {
      name = argv[0];
    }
    std::cerr << "Usage: " << name << " xyxYXY xxxYYYY\n";
  }

  Word u(argv[1]);
  Word v(argv[2]);

  auto canonical = GetCanonicalPair(u, v);

  std::cout << canonical.first << " " << canonical.second << std::endl;

  return 0;
}