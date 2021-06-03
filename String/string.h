//
// Created by Nikita Mastinen.
//

#include <cstring>
#include <iostream>

class String {
private:
  size_t real_size = 0;
  size_t buffer_size = 0;
  char* string_array= nullptr;
public:
  ~String() {
    delete[] string_array;
  }

  String() = default;

  String(const char* other) :
      real_size(strlen(other)),
      buffer_size(real_size),
      string_array(new char[buffer_size]) {
    memcpy(string_array, other, buffer_size);
  }

  String(size_t length, char character) :
      real_size(length),
      buffer_size(length),
      string_array(new char[buffer_size]) {
    memset(string_array, character, length);
  }

  String(char other) : String(1, other) {}

  String(const String& other) :
      real_size(other.real_size),
      buffer_size(other.buffer_size),
      string_array(new char[other.buffer_size]) {
    memcpy(string_array, other.string_array, buffer_size);
  }

  void swap(String& other) {
    std::swap(string_array, other.string_array);
    std::swap(real_size, other.real_size);
    std::swap(buffer_size, other.buffer_size);
  }

  String& operator=(const String& other) {
    if (other == (*this)) {
      return *this;
    }
    String copy = other;
    swap(copy);
    return *this;
  }

  char& operator[](size_t position) {
    return *(string_array + position);
  }

  char operator[](size_t position) const {
    return *(string_array + position);
  }

  size_t length() const {
    return real_size;
  }

  char& front() {
    return *(string_array);
  }

  char front() const {
    return *(string_array);
  }

  char& back() {
    return *(string_array + (real_size - 1));
  }

  char back() const {
    return *(string_array + (real_size - 1));
  }

  String& operator+=(char character) {
    if (buffer_size <= real_size) {
      buffer_size = (buffer_size == 0 ? 1 : 2 * buffer_size);
      char* copy = new char[buffer_size];
      memcpy(copy, string_array, real_size);
      delete[] string_array;
      string_array = copy;
    }
    *(string_array + real_size) = character;
    real_size++;
    return *this;
  }

  String& operator+=(const String& other) {
    size_t length = other.real_size;
    for (size_t i = 0; i < length; ++i) *this += other[i];
    return *this;
  }

  void push_back(char character) {
    *this += character;
  }

  void pop_back() {
    --real_size;
  }

  bool operator==(const String& other) const {
    if (real_size != other.real_size) {
      return false;
    }
    for (size_t i = 0; i < real_size; ++i) {
      if (string_array[i] != other.string_array[i]) {
        return false;
      }
    }
    return true;
  }

  bool empty() const {
    return real_size == 0;
  }

  void clear() {
    delete[] string_array;
    string_array = new char[0];
    real_size = 0;
    buffer_size = 0;
  }

  String substr(size_t start, size_t count) const {
    char* other = new char[count];
    memcpy(other, string_array + start, count);
    String copy;
    copy.real_size = count;
    copy.buffer_size = count;
    copy.string_array= other;
    return copy;
  }

  size_t find(const String& substring) const {
    for (size_t i = 0; i + substring.length() <= real_size; ++i) {
      if (substr(i, substring.length()) == substring) {
        return i;
      }
    }
    return real_size;
  }

  size_t rfind(const String& substring) const {
    if (real_size < substring.length()) return real_size;
    for (size_t i = real_size - substring.length(); ; --i) {
      if (substr(i, substring.length()) == substring) {
        return i;
      }
      if (i == 0) break;
    }
    return real_size;
  }
};

String operator+(const String& first_string, const String& second_string) {
  String sum = first_string;
  sum += second_string;
  return sum;
}

std::ostream& operator<<(std::ostream& out, const String& output) {
  for (size_t i = 0; i < output.length(); ++i) {
    out << output[i];
  }
  return out;
}

std::istream& operator>>(std::istream& in, String& input) {
  char character;
  input.clear();
  while (in >> std::noskipws >> character && character != ' ' && character != '\n') {
    input += character;
  }
  return in;
}
