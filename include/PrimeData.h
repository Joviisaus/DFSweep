#ifndef __PRIMEDATA_H__
#define __PRIMEDATA_H__
#include <vector>
struct PrimeData {
  int id;
  std::vector<double> params;
  int rank;
  double residual;
  bool isPlane;
};
#endif // !__PRIMEDATA_H__
