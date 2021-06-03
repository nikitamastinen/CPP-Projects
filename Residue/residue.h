//
// Created by Nikita Mastinen on 10.12.2020.
//

#include <cmath>
#include <iostream>
#include <vector>
#include <random>

using std::cin;
using std::cout;
using std::endl;
using std::vector;

template<unsigned long long N>
struct sqrt_calculator {
  static const unsigned root = static_cast<unsigned>(sqrt(N));
};

template<unsigned N>
const unsigned sqrt_calculator_v = sqrt_calculator<N>::root;

template<unsigned N, unsigned M>
struct is_prime_calculator {
  static const bool value = (N % M != 0 &&
                             is_prime_calculator<N, M - 1>::value);
};

template<unsigned N>
struct is_prime_calculator<N, 1> {
  static const bool value = true;
};

template<unsigned N>
struct is_prime {
  static const bool value = is_prime_calculator<N,
      sqrt_calculator_v<N>>::value;
};

template<>
struct is_prime<1> {
  static const bool value = false;
};

template<>
struct is_prime<0> {
  static const bool value = false;
};

template<unsigned N>
const bool is_prime_v = is_prime<N>::value;

template<unsigned N, unsigned M>
struct prime_divisor {
  static const unsigned long long value = (N % M == 0 && is_prime_v<M> ?
                                           M : prime_divisor<N, M - 1>::value);
};

template<unsigned N>
struct prime_divisor<N, 1> {
  static const unsigned long long value = N;
};

template<unsigned N>
struct prime_divisor<N, 2> {
  static const unsigned long long value = N;
};

template<unsigned N, unsigned M>
struct is_single_degree {
  static const unsigned long long value =
      M % 2 != 0 && (N % M == 0 ? is_single_degree<N / M, M>::value : (N == 1));
};

template<unsigned M>
struct is_single_degree<1, M> {
  static const unsigned long long value = true;
};

template<unsigned M>
struct is_single_degree<0, M> {
  static const unsigned long long value = false;
};

template<unsigned N>
struct has_primitive_root {
  static const bool value =
      is_single_degree<N / (2 - N % 2), prime_divisor<N / (2 - N % 2), sqrt_calculator_v<N / (2 - N % 2)>>::value>::value;
};

template<>
struct has_primitive_root<2> {
  static const bool value = true;
};

template<>
struct has_primitive_root<4> {
  static const bool value = true;
};

template<>
struct has_primitive_root<8> {
  static const bool value = false;
};

template<>
struct has_primitive_root<1> {
  static const bool value = false;
};

template<unsigned N>
const bool has_primitive_root_v = has_primitive_root<N>::value;

template<bool N>
struct throw_compilation_error {};

template<>
struct throw_compilation_error<false> {
  static const bool error = false;
};

int gcd(int x, int y) {
  return (y == 0 ? x : gcd(y, x % y));
}

unsigned euler_function(unsigned n) {
  int result = static_cast<int>(n);
  for (int i = 2; i * i <= n; i++) {
    if (n % i == 0) {
      while (n % i == 0)
        n /= i;
      result -= result / i;
    }
  }
  if (n > 1) {
    result -= result / n;
  }
  return result;
}


template<unsigned N>
class Residue {
  unsigned long long x = 0;
  static const unsigned euler_function_v ;
public:
  Residue() = default;

  explicit Residue(int x): x((x % static_cast<int>(N) + N) % N) {};

  explicit operator int() const {
    return x;
  }

  Residue& operator+=(Residue other) {
    x += other.x;
    if (x > N) x -= N;
    return *this;
  }

  Residue& operator*=(Residue other) {
    x *= other.x;
    x %= static_cast<long long>(N);
    return *this;
  }

  Residue& operator-=(Residue other) {
    x = x + N - other.x;
    if (x > N) x -= N;
    return *this;
  }

  Residue& operator/=(Residue other) {
    *this *= other.getInverse();
    return *this;
  }

  Residue pow(unsigned long long k) const {
    long long res = 1, t = x;
    while (k > 0) {
      if (k & 1) {
        res = res * t % static_cast<long long>(N);
      }
      t = t * t % static_cast<long long>(N);
      k /= 2;
    }
    return Residue(res);
  }

  Residue getInverse() const {
    if (throw_compilation_error<!is_prime_v<N>>::error) return Residue(0);
    return this->pow(N - 2);
  }

  unsigned int order() const {
    unsigned int update = 0;
    for (int i = 1; i * i <= euler_function_v; i++) {
      if (euler_function_v % i != 0) continue;
      if (pow(i).x == 1) {
        return i;
      }
      if (pow(euler_function_v / i).x == 1) {
        update = euler_function_v / i;
      }
    }
    return update;
  }

  static Residue getPrimitiveRoot() {
    for (int i = 2; i < static_cast<int>(N); i++) {
      if (gcd(euler_function_v, i) == 1 &&
          Residue(i).pow(euler_function_v / 2).x != 1) {
        bool is_primitive = true;
        for (int j = 2; j * j <= euler_function_v; j++) {
          if (euler_function_v % j == 0 && (Residue(i).pow(j).x == 1 || Residue(i).pow(euler_function_v / j).x == 1)) {
            is_primitive = false;
            break;
          }
        }
        if (is_primitive) return Residue(i);
      }
    }
    return Residue(0);
  }
};

template<unsigned N>
const unsigned Residue<N>::euler_function_v = euler_function(N);

template<unsigned N>
Residue<N> operator+(Residue<N> x, Residue<N> y) {
  x += y;
  return x;
}

template<unsigned N>
Residue<N> operator-(Residue<N> x, Residue<N> y) {
  x -= y;
  return x;
}

template<unsigned N>
Residue<N> operator*(Residue<N> x, Residue<N> y) {
  x *= y;
  return x;
}

template<unsigned N>
Residue<N> operator/(Residue<N> x, Residue<N> y) {
  x /= y;
  return x;
}