//
// Created by Nikita Mastinen on 23.10.2020.
//

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>
#include <complex>
#include <vector>
#include <algorithm>

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

class BigInteger;

BigInteger operator+(const BigInteger&, const BigInteger&);

BigInteger operator*(const BigInteger&, const BigInteger&);

BigInteger operator/(const BigInteger&, const BigInteger&);

class BigInteger {
  std::vector<short> integer = {};
  char sign = 0;

  void trim() {
    while(!integer.empty() && integer.back() == 0) {
      integer.pop_back();
    }
    if (integer.empty()) {
      sign = 0;
    }
  }

  void normalize() {
    for (size_t i = 0; i < integer.size(); ++i) {
      if (integer[i] >= 10) {
        if (i + 1 == integer.size()) {
          integer.push_back(integer[i] / 10);
        } else {
          integer[i + 1] += integer[i] / 10;
        }
        integer[i] %= 10;
      } else if (integer[i] < 0) {
        integer[i + 1] += integer[i] / 10;
        integer[i] %= 10;
        if (integer[i] != 0) {
          integer[i] += 10;
          integer[i + 1]--;
        }
      }
    }
  }

  void fast_normalize() {
    integer.push_back(0);
    for (size_t i = 0; i < integer.size(); ++i) {
      if (integer[i] >= 10) {
        ++integer[i + 1];
        integer[i] -= 10;
      } else if (integer[i] < 0) {
        --integer[i + 1];
        integer[i] += 10;
      }
    }
  }

  bool absolute_lower(const BigInteger& other) const {
    if (integer.size() != other.integer.size())
      return integer.size() < other.integer.size();
    for (size_t i = 0; i < integer.size(); ++i) {
      if (integer[integer.size() - 1 - i] < other.integer[integer.size() - 1 - i]) return true;
      if (integer[integer.size() - 1 - i] > other.integer[integer.size() - 1 - i]) return false;
    }
    return false;
  }

  BigInteger& add(const BigInteger& other, bool is_subtract) {
    if (other.integer.size() > integer.size()) {
      integer.resize(other.integer.size(), 0);
    }
    if ((!is_subtract ^ (sign != other.sign)) || sign == 0 || other.sign == 0) {
      if (sign == 0) sign = (is_subtract ? -other.sign : other.sign);
      for (size_t i = 0; i < other.integer.size(); ++i) {
        integer[i] += other.integer[i];
      }
    } else if (absolute_lower(other)) {
      if (!is_subtract) sign = other.sign;
      if (is_subtract) sign = -sign;
      for (size_t i = 0; i < other.integer.size(); ++i) {
        integer[i] = other.integer[i] - integer[i];
      }
    } else {
      if (!is_subtract) sign = -other.sign;
      for (size_t i = 0; i < other.integer.size(); ++i) {
        integer[i] -= other.integer[i];
      }
    }
    fast_normalize();
    trim();
    return *this;
  }

  static void push_front(std::vector<short>& v, short value) {
    v.push_back(0);
    for (size_t i = v.size() - 1; i > 0; --i) {
      v[i] = v[i - 1];
    }
    v[0] = value;
  }

public:
  BigInteger() = default;

  BigInteger(long long other) {
    if (other < 0) {
      sign = -1;
      other = -other;
    } else if (other == 0) {
      sign = 0;
    } else {
      sign = 1;
    }
    while (other > 0) {
      integer.push_back(other % 10);
      other /= 10;
    }
  }

  BigInteger(const std::string& other) {
    long long store = 0;
    sign = 1;
    size_t last = other.size();
    for (size_t i = 0; i < other.size(); ++i) {
      char character = other[other.size() - 1 - i];
      size_t is_signed = 0;
      if (character == '-') {
        sign = -1;
        is_signed = 1;
      }
      if (character == '+') {
        sign = 1;
        is_signed = 1;
      }
      if (other.size() - 1 - i == 0 || other.size() - i == last) {
        integer.push_back(0);
        for (size_t j = (last > 1 ? last - 1 : is_signed); j < last; ++j) {
          integer.back() = integer.back() * 10 + other[j] - '0';
        }
        --last;
      }
    }
    integer.push_back(store);
    trim();
    if (integer.empty()) {
      sign = 0;
    }
  }

  BigInteger(const BigInteger& other) = default;

  std::string toString() const {
    if (sign == 0) {
      return "0";
    }
    std::string result = (sign == -1 ? "-" : "");
    for (size_t i = 0; i < integer.size(); ++i) {
      if (i == 0) {
        result += std::to_string(integer[integer.size() - i - 1]);
      } else {
        int number_of_nulls =
            1 - static_cast<int>(std::to_string(integer[integer.size() - i - 1]).size());
        result += std::string(number_of_nulls, '0');
        result += std::to_string(integer[integer.size() - i - 1]);
      }
    }
    return result;
  }

  explicit operator bool() const {
    return sign != 0;
  }

  BigInteger& operator=(const BigInteger& other) = default;

  bool operator==(const BigInteger& other) const {
    return (sign == other.sign && integer == other.integer);
  }

  bool operator!=(const BigInteger& other) const {
    return !(*this == other);
  }

  bool operator<(const BigInteger& other) const {
    if (sign == other.sign) {
      return (sign == 1 ? absolute_lower(other) : !absolute_lower(other) && *this != other);
    }
    return sign < other.sign;
  }

  bool operator>=(const BigInteger& other) const {
    return !(*this < other);
  }

  bool operator<=(const BigInteger& other) const {
    return other >= *this;
  }

  bool operator>(const BigInteger& other) const {
    return other < *this;
  }

  BigInteger operator-() const {
    BigInteger negation = *this;
    negation.sign *= -1;
    return negation;
  }

  BigInteger& operator+=(const BigInteger& other) {
    return add(other, 0);
  }

  BigInteger& operator++() {
    *this += 1;
    return *this;
  }

  BigInteger operator++(int) {
    const BigInteger copy = *this;
    *this += 1;
    return copy;
  }

  BigInteger& operator-=(const BigInteger& other) {
    return add(other, true);
  }

  BigInteger& operator*=(const BigInteger& other) {
    BigInteger copy = *this;
    integer.clear();
    integer.resize(copy.integer.size() + other.integer.size());
    sign = copy.sign * other.sign;
    for (size_t i = 0; i < copy.integer.size(); ++i) {
      for (size_t j = 0; j < other.integer.size(); ++j) {
        integer[i + j] += copy.integer[i] * other.integer[j];
      }
    }
    normalize();
    trim();
    return *this;
  }

  int mod2() const {
    return integer.front() % 2;
  }

  int length() const {
    return integer.size();
  }

  BigInteger& div2() {
    for (size_t i = 0; i < integer.size(); ++i) {
      if (!(integer[i] & 1)) {
        integer[i] /= 2;
      } else {
        integer[i] /= 2;
        if (i) integer[i - 1] += 5;
      }
    }
    normalize();
    trim();
    return *this;
  }

  BigInteger& operator/=(BigInteger other) {
    if (other == 0) return *this;
    if (absolute_lower(other)) {
      *this = 0;
      return *this;
    }
    char new_sign = sign * other.sign;
    BigInteger copy = *this;
    copy.sign = ::abs(copy.sign);
    other.sign = ::abs(other.sign);
    *this = 0;
    BigInteger store = 0;
    store.sign = 1;
    store.integer.resize(other.integer.size());
    for (size_t i = other.integer.size(); i > 0; i--) {
      store.integer[i - 1] = (copy.integer[i - 1 - other.integer.size() + copy.integer.size()]);
    }
    int index = static_cast<int>(copy.integer.size() - other.integer.size());
    while (index >= 0) {
      bool full = false;
      while (index >= 1 && store < other) {
        --index;
        if (!(store.integer.empty() && copy.integer[index] == 0)) {
          push_front(store.integer, copy.integer[index]);
        }
        if (full) integer.push_back(0);
        full = true;
      }
      short cnt = 0;
      while (store >= other) {
        store -= other;
        cnt++;
      }
      store.sign = 1;
      integer.push_back(cnt);
      if (index <= 0) break;
    }
    std::reverse(integer.begin(), integer.end());
    sign = new_sign;
    trim();
    return *this;
  }

  BigInteger& operator%=(const BigInteger& other) {
    *this -= *this / other * other;
    return *this;
  }

  BigInteger abs() {
    BigInteger copy = *this;
    if (copy.sign < 0) {
      copy.sign = 1;
    }
    return copy;
  }
};

BigInteger operator*(const BigInteger& first_number, const BigInteger& second_number) {
  BigInteger result = first_number;
  result *= second_number;
  return result;
}

BigInteger operator+(const BigInteger& first_number, const BigInteger& second_number) {
  BigInteger result = first_number;
  result += second_number;
  return result;
}

BigInteger operator-(const BigInteger& first_number, const BigInteger& second_number) {
  BigInteger result = first_number;
  result -= second_number;
  return result;
}

BigInteger operator%(const BigInteger& first_number, const BigInteger& second_number) {
  BigInteger result = first_number;
  result %= second_number;
  return result;
}

BigInteger operator/(const BigInteger& first_number, const BigInteger& second_number) {
  BigInteger result = first_number;
  result /= second_number;
  return result;
}

BigInteger operator ""_bi(const unsigned long long other) {
  return BigInteger(static_cast<long long>(other));
}

std::istream& operator>>(std::istream& in, BigInteger& input) {
  std::string s;
  in >> s;
  input = s;
  return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& output) {
  out << output.toString();
  return out;
}

BigInteger gcd(BigInteger x, BigInteger y) {
  x = x.abs();
  y = y.abs();
  BigInteger res = 1;
  while (x != 0 && y != 0) {
    if (x.mod2() == 0 && y.mod2() == 1) x.div2();
    if (x.mod2() == 1 && y.mod2() == 0) y.div2();
    if (x.mod2() == 0 && y.mod2() == 0) x.div2(), y.div2(), res *= 2;
    if (x.mod2() == 1 && y.mod2() == 1) (x > y ? x = (x - y).div2() : y = (y - x).div2());
  }
  return (x == 0 ? y * res : x * res);
}

class Rational {
  BigInteger numerator = {};
  BigInteger denominator = {};

  void normalize() {
    if (denominator != 1) {
      BigInteger g = gcd(denominator, numerator);
      numerator /= g;
      denominator /= g;
    }
    if (denominator < 0) {
      denominator = -denominator;
      numerator = -numerator;
    }
  }
public:
  Rational(): numerator(0), denominator(1) {};

  Rational(const BigInteger& other): numerator(other), denominator(1) {};

  Rational(const long long other): numerator(other), denominator(1) {};

  Rational(const Rational& other) = default;;

  std::string toString() const {
    std::string output_numerator = numerator.toString();
    std::string output_denominator = denominator.toString();
    if (denominator != 1) {
      return output_numerator + "/" + output_denominator;
    } else {
      return output_numerator;
    }
  }

  Rational& operator=(const Rational& other) = default;

  bool operator==(const Rational& other) const {
    return numerator == other.numerator && denominator == other.denominator;
  }

  bool operator!=(const Rational& other) const {
    return !(*this == other);
  }

  bool operator<(const Rational& other) const {
    return numerator * other.denominator < denominator * other.numerator;
  }

  bool operator<=(const Rational& other) const {
    return numerator * other.denominator <= denominator * other.numerator;
  }

  bool operator>(const Rational& other) const {
    return numerator * other.denominator > denominator * other.numerator;
  }

  bool operator>=(const Rational& other) const {
    return numerator * other.denominator >= denominator * other.numerator;
  }

  Rational operator-() const {
    Rational negation = *this;
    negation.numerator = -negation.numerator;
    return negation;
  }

  Rational& operator+=(const Rational& other) {
    numerator = numerator * other.denominator + other.numerator * denominator;
    denominator *= other.denominator;
    normalize();
    return *this;
  }

  Rational& operator-=(const Rational& other) {
    numerator = numerator * other.denominator - other.numerator * denominator;
    denominator *= other.denominator;
    normalize();
    return *this;
  }

  Rational& operator*=(const Rational& other) {
    numerator *= other.numerator;
    denominator *= other.denominator;
    normalize();
    return *this;
  }

  Rational& operator/=(const Rational& other) {
    numerator *= other.denominator;
    denominator *= other.numerator;
    normalize();
    return *this;
  }

  std::string asDecimal(const size_t precision) const {
    Rational copy = *this;
    for (size_t i = 0; i < precision; ++i) {
      copy *= 10;
    }
    std::string current = (copy.numerator / copy.denominator).toString();
    std::string result;
    int is_negative = 0;
    if (current[0] == '-') {
      result += '-';
      is_negative = 2;
    }
    while (result.size() + current.size() - is_negative <= precision) {
      result += '0';
    }
    for (size_t i = (is_negative == 2 ? 1 : 0); i < current.size(); ++i) {
      result += current[i];
    }
    result.insert(result.begin() + result.size() - precision, '.');
    return result;
  }

  explicit operator double() const {
    std::stringstream copy;
    copy << asDecimal(100);
    double result;
    copy >> result;
    return result;
  }
};

Rational operator*(const Rational& first_number, const Rational& second_number) {
  Rational result = first_number;
  result *= second_number;
  return result;
}

Rational operator+(const Rational& first_number, const Rational& second_number) {
  Rational result = first_number;
  result += second_number;
  return result;
}

Rational operator-(const Rational& first_number, const Rational& second_number) {
  Rational result = first_number;
  result -= second_number;
  return result;
}

Rational operator/(const Rational& first_number, const Rational& second_number) {
  Rational result = first_number;
  result /= second_number;
  return result;
}

std::ostream& operator<<(std::ostream& out, const Rational& output) {
  out << output.toString();
  return out;
}


std::istream& operator>>(std::istream& in, Rational& input) {
  std::string s;
  in >> s;
  input = Rational(BigInteger(s));
  return in;
}

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
      is_single_degree<N / (2 - N % 2), prime_divisor<N / (2 - N % 2),
          sqrt_calculator_v<N / (2 - N % 2)>>::value>::value;
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
  for (int i = 2; i * i <= n; ++i) {
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

  Residue(int x): x((x % static_cast<int>(N) + N) % static_cast<int>(N)) {};

  explicit operator int() const {
    return x;
  }

  Residue& operator+=(Residue other) {
    x += other.x;
    if (x >= N) x -= N;
    return *this;
  }

  Residue& operator*=(Residue other) {
    x *= other.x;
    x %= static_cast<long long>(N);
    return *this;
  }

  bool operator==(Residue other) const {
    return x == other.x;
  }

  bool operator!=(Residue other) const {
    return x != other.x;
  }

  Residue& operator-=(Residue other) {
    x = x + N - other.x;
    if (x >= N) x -= N;
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
    for (int i = 1; i * i <= euler_function_v; ++i) {
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
    for (int i = 2; i < static_cast<int>(N); ++i) {
      if (gcd(euler_function_v, i) == 1 &&
          Residue(i).pow(euler_function_v / 2).x != 1) {
        bool is_primitive = true;
        for (int j = 2; j * j <= euler_function_v; ++j) {
          if (euler_function_v % j == 0 && (Residue(i).pow(j).x == 1 ||
                                            Residue(i).pow(euler_function_v / j).x == 1)) {
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
std::istream& operator>>(std::istream& in, Residue<N>& input) {
  int x;
  in >> x;
  input = Residue<N>(x);
  return in;
}

template<unsigned N>
std::ostream& operator<<(std::ostream& out, const Residue<N>& output) {
  out << static_cast<int>(output) << ",";
  return out;
}

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



template<unsigned N, unsigned M, typename Field>
class Matrix {
  vector<vector<Field>> matrix = vector<vector<Field>>(N, vector<Field>(M));
public:
  Matrix() = default;

  explicit Matrix(const vector<vector<Field>>& matrix): matrix(matrix) {};

  Matrix(const std::initializer_list<std::initializer_list<Field>>& list) {
    size_t index = 0;
    for (auto i : list) {
      size_t jndex = 0;
      for (auto j : i) {
        matrix[index][jndex] = Field(j);
        ++jndex;
      }
      ++index;
    }
  };

  Matrix& operator+=(const Matrix&);

  Matrix& operator-=(const Matrix&);

  Matrix& operator*=(const Field&);

  Matrix& operator*=(const Matrix&);

  template<unsigned N_new, unsigned M_new, unsigned K_new, typename Field_new>
  friend Matrix<N_new, K_new, Field_new> operator*
      (const Matrix<N_new, M_new, Field_new> &x,
       const Matrix<M_new, K_new, Field_new> &y);

  vector<Field> getRow(unsigned position) const;

  vector<Field> getColumn(unsigned position) const;

  unsigned rank() const;

  vector<Field>& operator[](unsigned position) {
    return matrix[position];
  }

  vector<Field>& operator[](unsigned position) const {
    return matrix[position];
  }

  Matrix<M, N, Field> transposed() const {
    vector<vector<Field>> result(M, vector<Field>(N));
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        result[j][i] = matrix[i][j];
      }
    }
    return Matrix<M, N, Field>(result);
  }

  bool operator==(const Matrix& other) const {
    return matrix == other.matrix;
  }

  bool operator!=(const Matrix& other) const {
    return matrix != other.matrix;
  }

  Field det() const;

  Field trace() const;

  void invert();

  Matrix inverted() const {
    auto copy = *this;
    copy.invert();
    return copy;
  }

  void swap(size_t i, size_t j) {
    std::swap(matrix[i], matrix[j]);
  }
};

template<unsigned N, unsigned M, typename Field = Rational>
void add(vector<vector<Field>>& x,
         const vector<vector<Field>>& y, bool is_subtraction) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      if (is_subtraction) x[i][j] -= y[i][j];
      else x[i][j] += y[i][j];
    }
  }
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator+=
    (const Matrix<N, M, Field>& other) {
  add<N, M, Field>(matrix, other.matrix, false);
  return *this;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator-=
    (const Matrix<N, M, Field>& other) {
  add<N, M, Field>(matrix, other.matrix, true);
  return *this;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=
    (const Field &value) {
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      matrix[i][j] *= value;
    }
  }
  return *this;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator*=
    (const Matrix& other) {
  static_assert(N == M);
  *this = (*this) * other;
  return *this;
}


template<unsigned N, unsigned M, typename Field>
vector<Field> Matrix<N, M, Field>::getRow(unsigned position) const {
  vector<Field> result;
  for (size_t i = 0; i < M; ++i) {
    result.push_back(matrix[position][i]);
  }
  return result;
}

template<unsigned N, unsigned M, typename Field>
vector<Field> Matrix<N, M, Field>::getColumn(unsigned position) const {
  vector<Field> result;
  for (size_t i = 0; i < N; ++i) {
    result.push_back(matrix[i][position]);
  }
  return result;
}

template<unsigned N, unsigned M, unsigned K, typename Field>
Matrix<N, K, Field> operator*(const Matrix<N, M, Field> &x, const Matrix<M, K, Field> &y) {
  Matrix<N, K, Field> result;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < K; ++j) {
      for (size_t k = 0; k < M; ++k) {
        result.matrix[i][j] += x.matrix[i][k] * y.matrix[k][j];
      }
    }
  }
  return result;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator*(const Field& value, const Matrix<N, M, Field> &y) {
  Matrix<N, M, Field> result = y;
  result *= value;
  return result;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> step(const Matrix<N, M, Field> &x) {
  cerr << N << endl;
  Matrix<N, M, Field> a = x;
  vector<bool> used(M);
  for (size_t k = 0; k < M; ++k) {
    size_t position = N + 1;
    for (size_t i = 0; i < N; ++i) {
      if (a[i][k] == Field(0) || used[i]) continue;
      if (position == N + 1) {
        position = i;
        used[i] = true;
      } else {
        Field value = a[i][k] / a[position][k];
        for (size_t j = k; j < M; ++j) {
          a[i][j] -= a[position][j] * value;
        }
      }
    }
  }
  return a;
}

template<unsigned N, unsigned M, typename Field>
unsigned Matrix<N, M, Field>::rank() const {
  unsigned result = 0;
  Matrix<N, M, Field> a = step(*this);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      if (a[i][j] != Field(0)) {
        ++result;
        break;
      }
    }
  }
  return result;
}

template<unsigned int N, unsigned int M, typename Field>
Field Matrix<N, M, Field>::det() const {
  static_assert(N == M);
  Field result = Field(1);
  Matrix<N, M, Field> a = step(*this);
  for (size_t i = 0; i < N; ++i) {
    result *= a[i][i];
  }
  return result;
}

template<unsigned int N, unsigned int M, typename Field>
Field Matrix<N, M, Field>::trace() const {
  static_assert(N == M);
  Field result = 0;
  for (size_t i = 0; i < N; ++i) {
    result += matrix[i][i];
  }
  return result;
}

template<unsigned int N, unsigned int M, typename Field>
void Matrix<N, M, Field>::invert() {
  static_assert(N == M);
  Matrix<N, N + N, Field> result;
  for (size_t i = 0; i < N; ++i) {
    result[i][i + N] = Field(1);
  }
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      result[i][j] = matrix[i][j];
    }
  }
  result = step(result);
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (result[i][j] != Field(0)) {
        result.swap(i, j);
        break;
      }
    }
  }
  for (size_t i = 0; i < N; ++i) {
    if (result[i][i] != Field(1)) {
      Field value = result[i][i];
      for (size_t j = i; j < N + N; ++j) {
        result[i][j] /= value;
      }
    }
  }
  for (size_t i = 0; i < N; ++i) {
    result.swap(i, N - i - 1);
  }
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = i + 1; j < N; ++j) {
      if (result[i][j] != Field(0)) {
        Field value = result[i][j];
        for (size_t k = j; k < N + N; ++k) {
          result[i][k] -= result[j][k] * value;
        }
      }
    }
  }
  for (size_t i = 0; i < N; ++i) {
    result.swap(i, N - i - 1);
  }
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      matrix[i][j] = result[i][N + j];
    }
  }
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator+
    (const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
  Matrix<N, M, Field> copy = x;
  copy += y;
  return copy;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator-
    (const Matrix<N, M, Field>& x, const Matrix<N, M, Field>& y) {
  Matrix<N, M, Field> copy = x;
  copy -= y;
  return copy;
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator*
    (const Matrix<N, M, Field>& x, const Field& value) {
  Matrix<N, M, Field> copy = x;
  copy *= value;
  return copy;
}

template<unsigned N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;