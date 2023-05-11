#pragma once

#include <compare>
#include <concepts>
#include <functional>
#include <iosfwd>
#include <limits>
#include <string>
#include <utility>
#include <vector>

struct big_integer {
  big_integer();
  big_integer(const big_integer& other);

  big_integer(unsigned long long a);
  big_integer(long long a);
  big_integer(unsigned int a);
  big_integer(int a);
  big_integer(unsigned long a);
  big_integer(long a);

  explicit big_integer(const std::string& str);
  ~big_integer();

  big_integer& operator=(const big_integer& other);

  big_integer& operator+=(const big_integer& rhs);
  big_integer& operator-=(const big_integer& rhs);
  big_integer& operator*=(const big_integer& rhs);
  big_integer& operator/=(const big_integer& rhs);
  big_integer& operator%=(const big_integer& rhs);

  big_integer& operator&=(const big_integer& rhs);
  big_integer& operator|=(const big_integer& rhs);
  big_integer& operator^=(const big_integer& rhs);

  big_integer& operator<<=(int rhs);
  big_integer& operator>>=(int rhs);

  big_integer operator+() const;
  big_integer operator-() const;
  big_integer operator~() const;

  big_integer& operator++();
  big_integer operator++(int);

  big_integer& operator--();
  big_integer operator--(int);

  friend std::strong_ordering operator<=>(const big_integer& a, const big_integer& b);
  friend bool operator==(const big_integer& a, const big_integer& b);
  friend bool operator!=(const big_integer& a, const big_integer& b);
  friend bool operator<(const big_integer& a, const big_integer& b);
  friend bool operator>(const big_integer& a, const big_integer& b);
  friend bool operator<=(const big_integer& a, const big_integer& b);
  friend bool operator>=(const big_integer& a, const big_integer& b);

  friend std::string to_string(const big_integer& a);

  friend big_integer operator*(const big_integer& a, const big_integer& b);
  friend big_integer operator/(const big_integer& a, const big_integer& b);
  friend big_integer operator%(const big_integer& a, const big_integer& b);

private:
  using digit = uint32_t;
  using double_digit = uint64_t;

  static const digit MAX_DIGIT;
  static const double_digit BASE;
  static const size_t DIGIT_SIZE;
  static const size_t EXP10;
  static const big_integer BASE10;

  big_integer(digit d, bool is_negative);
  big_integer(unsigned long long a, bool is_negative);

  std::vector<digit> digits_;
  bool is_negative_;

  big_integer& abs_add_digit(digit d);
  big_integer& abs_sub_digit(digit d);
  big_integer& add_digit(digit d);
  big_integer& sub_digit(digit d);
  big_integer& negate();
  big_integer& negate_if(bool cond);
  big_integer& to_abs();
  big_integer& to_zero();
  big_integer& abs_add_shifted(const big_integer& rhs, size_t shift = 0);
  big_integer& abs_sub_shifted(const big_integer& rhs, size_t shift = 0);
  big_integer& add_shifted(const big_integer& rhs, size_t shift = 0);
  big_integer& sub_shifted(const big_integer& rhs, size_t shift = 0);
  big_integer& bitwise(const big_integer& rhs, const std::function<digit(digit, digit)>& f);

  void swap(big_integer& other);
  void strip_zeros() noexcept;

  size_t size() const;

  std::strong_ordering abs_compare_shifted(const big_integer& rhs, size_t shift = 0) const;
  big_integer mul_digit(digit d) const;
  bool is_zero() const;
  size_t get_norm() const;
  std::pair<big_integer, big_integer> divrem(const big_integer& b) const;
  big_integer abs() const;
  void check_invariant() const;
  digit get_digit_value(size_t n) const;
};

big_integer operator+(const big_integer& a, const big_integer& b);
big_integer operator-(const big_integer& a, const big_integer& b);

big_integer operator&(const big_integer& a, const big_integer& b);
big_integer operator|(const big_integer& a, const big_integer& b);
big_integer operator^(const big_integer& a, const big_integer& b);

big_integer operator<<(const big_integer& a, int b);
big_integer operator>>(const big_integer& a, int b);

bool operator==(const big_integer& a, const big_integer& b);
bool operator!=(const big_integer& a, const big_integer& b);
bool operator<(const big_integer& a, const big_integer& b);
bool operator>(const big_integer& a, const big_integer& b);
bool operator<=(const big_integer& a, const big_integer& b);
bool operator>=(const big_integer& a, const big_integer& b);

std::string to_string(const big_integer& a);
std::ostream& operator<<(std::ostream& out, const big_integer& a);
