#include <bits/stdc++.h>
// #include <algorithm>
// using namespace std;

#define xstr(s) str(s)
#define str(s) #s
/*
#ifndef NELE
#define NELE 3
#endif
*/
#ifndef TYPE
#define TYPE double
#endif

#ifndef NSAMPLE
#define NSAMPLE 1
#endif

#ifndef OPTIM
#define OPTIM classic
#endif

#define STROPTIM xstr(OPTIM)
#define STRTYPE xstr(TYPE)

inline bool IsFileExist(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}

int main(int argc, char *argv[]) {
  int ierr = 0;
  long NELE = std::stol(argv[1]);
  std::vector<TYPE> A(NELE, 1);
  std::vector<TYPE> B(NELE, 2);
  std::vector<TYPE> C(NELE, 3);
  std::vector<TYPE> RES(NELE, 0);
  double Mem = 4 * NELE * sizeof(TYPE) / 1000000.; // MB
  auto t_before = std::chrono::steady_clock::now();
  for (long i = 0; i < NSAMPLE; i++) {
    std::transform(A.begin(), A.end(), B.begin(), RES.begin(),
                   std::multiplies<TYPE>());
    std::transform(C.begin(), C.end(), RES.begin(), RES.begin(),
                   std::plus<TYPE>());
  }
  auto t_after = std::chrono::steady_clock::now();
  auto d_delta = t_after - t_before;
  auto decod_time_ms = (float)d_delta.count() * 0.000001f;
  std::cout << STRTYPE << " " << STROPTIM << '\n';
  std::cout << "time: " << decod_time_ms / NSAMPLE << " ms" << std::endl;
  std::cout << "cumulative time: " << decod_time_ms << " ms over " << NSAMPLE
            << " Mem MB " << Mem << std::endl;
  TYPE sum = std::accumulate(RES.begin(), RES.end(), (double)0.0);
  std::cout << ((5 * (double)NELE) == sum) << " " << (5) * NELE << " " << sum
            << '\n';
  std::ofstream tabledb;
  auto fname = "benchfmadd.csv";
  auto header = "type,Total,Time,NSample,Mem,optim";
  auto isfile = IsFileExist(fname);
  if (!isfile) {
    tabledb.open(fname, std::ios::app);
    tabledb << header << "\n";
    tabledb.close();
  }
  tabledb.open(fname, std::ios::app);
  tabledb << STRTYPE << "," << decod_time_ms << "," << decod_time_ms / NSAMPLE
          << "," << NSAMPLE << "," << Mem << "," << STROPTIM << "\n";
  tabledb.close();
  return ierr;
}