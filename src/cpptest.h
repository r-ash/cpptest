#ifndef _CPPTEST_CPPTEST_H_
#define _CPPTEST_CPPTEST_H_

double add(double a, double b);

class doubler {
private:
  double number;
public:
  // doubler(double x) {
  //   number = x;
  // }
  doubler(double x) : number(x) {}
  double run() {
    number = number * 2;
    return number;
  }
};

#endif
