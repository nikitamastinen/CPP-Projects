//
// Created by Nikita Mastinen on 22.10.2020.
//

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <complex>

using complex = std::complex<double>;

size_t reverse_bits(size_t x, size_t n) {
  size_t result = 0;
  for (size_t i = 2, cnt = 0; i <= n; i *= 2, cnt++) {
    size_t bit = (bool)(n / i & x);
    result += (bit << cnt);
  }
  return result;
}

std::vector<complex> fft(std::vector<complex> polynom, bool reverse = false) {
  for (size_t i = 0; i < polynom.size(); i++) {
    if (i < reverse_bits(i, polynom.size())) {
      swap(polynom[i], polynom[reverse_bits(i, polynom.size())]);
    }
  }
  std::vector<complex> roots(polynom.size(), 1);
  for (size_t i = 0; i < polynom.size(); i++) {
    double angle = (reverse ? -1 : 1) * acos(-1) * 2 * i / polynom.size();
    roots[i] = complex(cos(angle), sin(angle));
  }
  for (size_t i = 2; i <= polynom.size(); i *= 2) {
    for (size_t j = 0; j < polynom.size(); j += i) {
      complex counter = 1;
      for (size_t k = j; k < j + i / 2; k++) {
        complex even = polynom[k] + polynom[k + i / 2] * counter;
        complex odd = polynom[k] - polynom[k + i / 2] * counter;
        polynom[k] = even;
        polynom[k + i / 2] = odd;
        counter = roots[polynom.size() / i * (k - j + 1)];
      }
    }
  }
  return polynom;
}

std::vector<int> multiply(const std::vector<int>& first_polynom,
                          const std::vector<int>& second_polynom) {
  size_t degree = 1;
  while (degree < std::max(first_polynom.size(), second_polynom.size())) degree *= 2;
  degree*=2;
  std::vector<complex> fft_first(degree), fft_second(degree), fft_result(degree);
  std::vector<int> result(degree);
  for (size_t i = 0; i < first_polynom.size(); ++i) fft_first[i] = first_polynom[i];
  for (size_t i = 0; i < second_polynom.size(); ++i) fft_second[i] = second_polynom[i];
  fft_first = fft(fft_first), fft_second = fft(fft_second);
  for (size_t i = 0; i < degree; i++) {
    fft_result[i] = fft_first[i] * fft_second[i];
  }
  fft_first.clear();
  fft_second.clear();
  fft_result = fft(fft_result, true);
  for (size_t i = 0; i < degree; i++) {
    result[i] = round(fft_result[i].real() / degree);
  }
  return result;
}

class BigInteger;

BigInteger operator+(const BigInteger& first_number, const BigInteger& second_number);

BigInteger operator*(const BigInteger& first_number, const BigInteger& second_number);

BigInteger operator/(const BigInteger& first_number, const BigInteger& second_number);

class BigInteger {
  static const int base = 10000;
  static const int number_of_digits = 4;
  std::vector<int> integer = {};
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
      if (integer[i] >= base) {
        if (i + 1 == integer.size()) {
          integer.push_back(integer[i] / base);
        } else {
          integer[i + 1] += integer[i] / base;
        }
        integer[i] %= base;
      } else if (integer[i] < 0) {
        integer[i + 1] += integer[i] / base;
        integer[i] %= base;
        if (integer[i] != 0) {
          integer[i] += base;
          integer[i + 1]--;
        }
      }
    }
  }

  bool absolute_lower(const BigInteger& other) const {
    if (integer.size() != other.integer.size())
      return integer.size() < other.integer.size();
    for (size_t i = 0; i < integer.size(); i++) {
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
      for (size_t i = 0; i < other.integer.size(); i++) {
        integer[i] += other.integer[i];
      }
    } else if (absolute_lower(other)) {
      if (!is_subtract) sign = other.sign;
      if (is_subtract) sign = -sign;
      for (size_t i = 0; i < other.integer.size(); i++) {
        integer[i] = other.integer[i] - integer[i];
      }
    } else {
      if (!is_subtract) sign = -other.sign;
      for (size_t i = 0; i < other.integer.size(); i++) {
        integer[i] -= other.integer[i];
      }
    }
    normalize();
    trim();
    return *this;
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
      integer.push_back(other % base);
      other /= base;
    }
  }

  BigInteger(const std::string& other) {
    long long store = 0;
    sign = 1;
    size_t last = other.size();
    for (size_t i = 0; i < other.size(); i++) {
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
      if (other.size() - 1 - i == 0 || other.size() - i - 1 + number_of_digits == last) {
        integer.push_back(0);
        for (size_t j = (last > number_of_digits ? last - number_of_digits : is_signed); j < last; j++) {
          integer.back() = integer.back() * 10 + other[j] - '0';
        }
        last -= number_of_digits;
      }
    }
    integer.push_back(store);
    trim();
    if (integer.empty()) {
      sign = 0;
    }
  }

  BigInteger(const BigInteger& other) : integer(other.integer), sign(other.sign) {}

  std::string toString() const {
    if (sign == 0) {
      return "0";
    }
    std::string result = (sign == -1 ? "-" : "");
    for (size_t i = 0; i < integer.size(); i++) {
      if (i == 0) {
        result += std::to_string(integer[integer.size() - i - 1]);
      } else {
        int number_of_nulls =
            number_of_digits - static_cast<int>(std::to_string(integer[integer.size() - i - 1]).size());
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
    char new_sign = sign * other.sign;
    integer = multiply(integer, other.integer);
    sign = new_sign;
    normalize();
    trim();
    return *this;
  }

  BigInteger& operator/=(BigInteger other) {
    if (other == 0) return *this;
    char new_sign = sign * other.sign;
    BigInteger copy = *this;
    copy.sign = ::abs(copy.sign);
    other.sign = ::abs(other.sign);
    *this = 0;
    while (other.absolute_lower(copy + 1)) {
      BigInteger store = other, factor = 1;
      while ((store * base).absolute_lower(copy)) {
        store *= base;
        factor *= base;
      }
      while (store.absolute_lower(copy + 1)) {
        copy -= store;
        *this += factor;
      }
    }
    sign = new_sign;
    normalize();
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

BigInteger gcd(const BigInteger& x, const BigInteger& y) {
  return (y == 0 ? x : gcd(y, x % y));
}

class Rational {
  BigInteger numerator = {};
  BigInteger denominator = {};

  void normalize() {
    BigInteger g = gcd(denominator, numerator).abs();
    numerator /= g;
    denominator /= g;
    if (denominator < 0) {
      denominator = -denominator;
      numerator = -numerator;
    }
  }
public:
  Rational(): numerator(0), denominator(1) {};

  Rational(const BigInteger& other): numerator(other), denominator(1) {};

  Rational(const long long other): numerator(other), denominator(1) {};

  Rational(const Rational& other):
      numerator(other.numerator),
      denominator(other.denominator) {
  };

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
    return other < *this;
  }

  bool operator>=(const Rational& other) const {
    return other <= *this;
  }

  Rational operator-() const {
    Rational negation = *this;
    negation.numerator = -negation.numerator;
    return negation;
  }

  Rational& operator+=(const Rational& other) {
    numerator = numerator * other.denominator + other.numerator * denominator;
    denominator = denominator * other.denominator;
    normalize();
    return *this;
  }

  Rational& operator-=(const Rational& other) {
    numerator = numerator * other.denominator - other.numerator * denominator;
    denominator = denominator * other.denominator;
    normalize();
    return *this;
  }

  Rational& operator*=(const Rational& other) {
    numerator = numerator * other.numerator;
    denominator = denominator * other.denominator;
    normalize();
    return *this;
  }

  Rational& operator/=(const Rational& other) {
    if (*this == other) {
      *this = 1;
      return *this;
    }
    numerator = numerator * other.denominator;
    denominator = denominator * other.numerator;
    normalize();
    return *this;
  }

  std::string asDecimal(const size_t precision) const {
    Rational copy = *this;
    for (size_t i = 0; i < precision; i++) {
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
    for (size_t i = (is_negative == 2 ? 1 : 0); i < current.size(); i++) {
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