#pragma once

#include <compare>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

struct big_integer {
  big_integer();
  big_integer(const big_integer& other);
  big_integer(int a);
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

    using digit = uint32_t;
    using double_digit = uint64_t;
      static const digit base;
    
  big_integer& negate();
  big_integer& negate_if(bool cond);
  //big_integer& inverse();
  size_t size() const;

  friend std::string to_string(const big_integer& a);


  friend big_integer operator*(const big_integer& a, const big_integer& b);
  friend big_integer operator/(const big_integer& a, const big_integer& b);
  friend big_integer operator%(const big_integer& a, const big_integer& b);

  private:
  static const size_t digit_size;
  static const size_t exp10;
  static const big_integer base10;

  std::vector<digit> digits_;
  bool is_negative_;

  void swap(big_integer& other);
  big_integer& bitwise(const big_integer& rhs, std::function<digit(digit, digit)> f);
  big_integer mul_digit(digit d) const;
  void strip_zeros();
  std::strong_ordering abs_compare_shifted(const big_integer& rhs, size_t shift = 0) const;
  bool is_zero() const;
  size_t get_norm() const;
  std::pair<big_integer, big_integer> divrem(const big_integer& b) const;
  big_integer& to_abs();
  big_integer abs() const;

  big_integer& abs_add_shifted(const big_integer& rhs, size_t shift = 0);
  big_integer& abs_sub_shifted(const big_integer& rhs, size_t shift = 0);
  big_integer& add_shifted(const big_integer& rhs, size_t shift = 0);
  big_integer& sub_shifted(const big_integer& rhs, size_t shift = 0);

  big_integer(digit d, bool is_negative);
};

big_integer operator+(const big_integer& a, const big_integer& b);
big_integer operator-(const big_integer& a, const big_integer& b);
/*
big_integer operator*(const big_integer& a, const big_integer& b);
big_integer operator/(const big_integer& a, const big_integer& b);
big_integer operator%(const big_integer& a, const big_integer& b);
*/

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
