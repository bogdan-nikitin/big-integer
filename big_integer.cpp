#include "big_integer.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <compare>
#include <cstddef>
#include <functional>
#include <limits>
#include <ostream>
#include <stack>
#include <stdexcept>
#include <string>

const size_t big_integer::digit_size = std::numeric_limits<big_integer::digit>::digits;
const big_integer::digit big_integer::base = std::numeric_limits<big_integer::digit>::max();
const size_t big_integer::exp10 = std::numeric_limits<digit>::digits10;
const big_integer big_integer::base10(static_cast<digit>(std::pow(10, big_integer::exp10)), false);

static std::strong_ordering less(bool cond) {
  return cond ? std::strong_ordering::less : std::strong_ordering::greater;
}

static std::strong_ordering reverse(std::strong_ordering order) {
  return 0 <=> order;
}

big_integer::big_integer() : is_negative_{false} {}

big_integer::big_integer(const big_integer& other) = default;

#define BIG_INTEGER_CTOR(T)                                                                                            \
  big_integer::big_integer(T a) { from_primitive(a); }

BIG_INTEGER_CTOR(short)
BIG_INTEGER_CTOR(unsigned short)
BIG_INTEGER_CTOR(unsigned)
BIG_INTEGER_CTOR(int)
BIG_INTEGER_CTOR(long)
BIG_INTEGER_CTOR(unsigned long)
BIG_INTEGER_CTOR(long long)
BIG_INTEGER_CTOR(unsigned long long)

#undef BIG_INTEGER_CTOR

big_integer::big_integer(const std::string& str) : is_negative_(str.starts_with("-") && str != "-0") {
  if (str == "-0") {
    return;
  }
  if (str.empty() || std::any_of(str.begin() + is_negative_, str.end(), [](char c) { return !std::isdigit(c); })) {
    throw std::invalid_argument("Invalid string");
  }
  size_t bound = exp10 + is_negative_;
  size_t i = str.size();
  big_integer cur_base(1);
  while (true) {
    bool in_bound = i >= bound;
    size_t begin = in_bound ? i - exp10 : is_negative_;
    size_t end = in_bound ? exp10 : i - is_negative_;
    abs_add_shifted(cur_base * big_integer(static_cast<digit>(std::stoul(str.substr(begin, end))), false));
    i = in_bound ? i - exp10 : is_negative_;
    if (i == is_negative_) {
      break;
    }
    cur_base *= base10;
  }
  check_invariant();
}

big_integer::~big_integer() = default;

big_integer& big_integer::operator=(const big_integer& other) = default;

big_integer& big_integer::operator+=(const big_integer& rhs) {
  return add_shifted(rhs);
}

big_integer& big_integer::operator-=(const big_integer& rhs) {
  return sub_shifted(rhs);
}

big_integer& big_integer::operator*=(const big_integer& rhs) {
  (*this * rhs).swap(*this);
  return *this;
}

big_integer& big_integer::operator/=(const big_integer& rhs) {
  (*this / rhs).swap(*this);
  return *this;
}

big_integer& big_integer::operator%=(const big_integer& rhs) {
  (*this % rhs).swap(*this);
  return *this;
}

big_integer& big_integer::operator&=(const big_integer& rhs) {
  return bitwise(rhs, std::bit_and());
}

big_integer& big_integer::operator|=(const big_integer& rhs) {
  return bitwise(rhs, std::bit_or());
}

big_integer& big_integer::operator^=(const big_integer& rhs) {
  return bitwise(rhs, std::bit_xor());
}

big_integer& big_integer::operator<<=(int rhs) {
  auto shift = static_cast<digit>(rhs);
  size_t start = shift / digit_size;
  size_t gain = start + (shift % digit_size != 0);
  size_t remain = shift - start * digit_size;
  size_t old_size = size();
  digits_.resize(old_size + gain);
  if (remain == 0) {
    std::shift_right(digits_.begin(), digits_.end(), static_cast<ptrdiff_t>(gain));
  } else {
    digit prev = 0;
    for (size_t i = old_size; i > 0; --i) {
      digits_[i - 1 + gain] = (prev << remain) | (digits_[i - 1] >> (digit_size - remain));
      prev = digits_[i - 1];
    }
    digits_[start] = digits_[0] << remain;
  }
  std::fill(digits_.begin(), digits_.begin() + static_cast<ptrdiff_t>(start), 0);
  strip_zeros();
  check_invariant();
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
  auto shift = static_cast<digit>(rhs);
  size_t loss = shift / digit_size;
  if (loss >= size()) {
    return to_zero();
  }
  size_t remain = shift - loss * digit_size;
  bool is_neg = is_negative_;
  if (is_neg) {
    ++*this;
  }
  if (remain == 0) {
    std::shift_left(digits_.begin(), digits_.end(), static_cast<ptrdiff_t>(loss));
  } else {
    size_t last = size() - loss - 1;
    for (size_t i = 0; i < last; ++i) {
      digits_[i] = (digits_[i + loss] >> remain) | (digits_[i + 1 + loss] << (digit_size - remain));
    }
    digits_[last] = digits_[last + loss] >> remain;
  }
  digits_.resize(size() - loss);
  strip_zeros();
  if (is_neg) {
    --*this;
  }
  check_invariant();
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  return big_integer(*this).negate();
}

big_integer big_integer::operator~() const {
  return --big_integer{*this}.negate();
}

big_integer& big_integer::operator++() {
  if (is_negative_) {
    bool borrow = true;
    for (size_t i = 0; i < size() && borrow; ++i) {
      borrow = digits_[i] == 0;
      digits_[i]--;
    }
    strip_zeros();
  } else {
    bool carry = true;
    for (size_t i = 0; i < size() && carry; ++i) {
      digits_[i]++;
      carry = digits_[i] == 0;
    }
    if (carry) {
      digits_.push_back(1);
    }
  }
  check_invariant();
  return *this;
}

big_integer big_integer::operator++(int) {
  big_integer result{*this};
  ++*this;
  return result;
}

big_integer& big_integer::operator--() {
  return (++negate()).negate();
}

big_integer big_integer::operator--(int) {
  big_integer result{*this};
  --*this;
  return result;
}

void big_integer::swap(big_integer& other) {
  std::swap(digits_, other.digits_);
  std::swap(is_negative_, other.is_negative_);
}

big_integer operator+(const big_integer& a, const big_integer& b) {
  return big_integer(a) += b;
}

big_integer operator-(const big_integer& a, const big_integer& b) {
  return big_integer(a) -= b;
}

big_integer operator*(const big_integer& a, const big_integer& b) {
  big_integer result;
  for (size_t i = 0; i < b.size(); ++i) {
    result.abs_add_shifted(a.mul_digit(b.digits_[i]), i);
  }
  result.is_negative_ = a.is_negative_ ^ b.is_negative_;
  result.check_invariant();
  return result;
}

big_integer operator/(const big_integer& a, const big_integer& b) {
  return a.divrem(b).first;
}

big_integer operator%(const big_integer& a, const big_integer& b) {
  return a.divrem(b).second;
}

big_integer operator&(const big_integer& a, const big_integer& b) {
  return big_integer(a) &= b;
}

big_integer operator|(const big_integer& a, const big_integer& b) {
  return big_integer(a) |= b;
}

big_integer operator^(const big_integer& a, const big_integer& b) {
  return big_integer(a) ^= b;
}

big_integer operator<<(const big_integer& a, int b) {
  return big_integer(a) <<= b;
}

big_integer operator>>(const big_integer& a, int b) {
  return big_integer(a) >>= b;
}

std::strong_ordering operator<=>(const big_integer& a, const big_integer& b) {
  if (a.is_negative_ != b.is_negative_) {
    return less(a.is_negative_);
  }
  std::strong_ordering order = a.abs_compare_shifted(b);
  return a.is_negative_ ? reverse(order) : order;
}

bool operator==(const big_integer& a, const big_integer& b) = default;
bool operator!=(const big_integer& a, const big_integer& b) = default;
bool operator<(const big_integer& a, const big_integer& b) = default;
bool operator>(const big_integer& a, const big_integer& b) = default;
bool operator<=(const big_integer& a, const big_integer& b) = default;
bool operator>=(const big_integer& a, const big_integer& b) = default;

big_integer& big_integer::negate() {
  is_negative_ = !is_zero() && !is_negative_;
  check_invariant();
  return *this;
}

big_integer& big_integer::negate_if(bool cond) {
  is_negative_ = cond && !is_zero() && !is_negative_;
  check_invariant();
  return *this;
}

big_integer& big_integer::bitwise(
    const big_integer& rhs, 
    const std::function<big_integer::digit(big_integer::digit, big_integer::digit)>& f) {
  if (size() < rhs.size()) {
    digits_.resize(rhs.size());
  }
  bool borrow = is_negative_, rhs_borrow = rhs.is_negative_;
  bool is_neg = f(is_negative_, rhs.is_negative_);
  size_t i = 0;
  for (; i < rhs.size(); ++i) {
    digit d = digits_[i] - borrow;
    digit rhs_d = rhs.digits_[i] - rhs_borrow;
    borrow = digits_[i] == 0 && borrow;
    rhs_borrow = rhs.digits_[i] == 0 && rhs_borrow;
    digits_[i] = f(is_negative_ ? ~d : d, rhs.is_negative_ ? ~rhs_d : rhs_d) ^ (is_neg ? big_integer::base : 0);
  }
  for (; i < size(); ++i) {
    digit d = digits_[i] - borrow;
    digits_[i] = f(is_negative_ ? ~d : d, (rhs.is_negative_ ? base : 0)) ^ (is_neg ? big_integer::base : 0);
  }
  is_negative_ = is_neg;
  strip_zeros();
  if (is_neg) {
    --*this;
  }
  check_invariant();
  return *this;
}

big_integer big_integer::mul_digit(digit d) const {
  if (d == 0 || is_zero()) {
    return {};
  }
  big_integer result;
  result.digits_.resize(size());
  bool carry = false;
  for (size_t i = 0; i < size() - 1; ++i) {
    double_digit product = static_cast<double_digit>(digits_[i]) * d;
    digit lo = product;
    digit& rd = result.digits_[i];
    digit sum = rd + product + carry;
    carry = rd > base - lo || (rd + lo == base && carry);
    rd = sum;
    result.digits_[i + 1] += product >> digit_size;
  }
  double_digit product = static_cast<double_digit>(digits_[size() - 1]) * d;
  digit lo = product;
  digit& rd = result.digits_[size() - 1];
  digit sum = rd + product + carry;
  carry = rd > base - lo || (rd + lo == base && carry);
  rd = sum;
  digit hi = product >> digit_size;
  if (hi != 0 || carry) {
    result.digits_.push_back(hi + carry);
    if (hi == base && carry) {
      result.digits_.push_back(1);
    }
  }
  result.check_invariant();
  return result;
}

size_t big_integer::size() const {
  return digits_.size();
}

void big_integer::strip_zeros() {
  size_t i = size();
  while (i > 0 && digits_[i - 1] == 0) {
    --i;
  }
  digits_.erase(digits_.begin() + static_cast<ptrdiff_t>(i), digits_.end());
  is_negative_ &= !is_zero();
  check_invariant();
}

std::strong_ordering big_integer::abs_compare_shifted(const big_integer& rhs, size_t shift) const {
  size_t rhs_size = rhs.size() + shift;
  if (size() != rhs_size) {
    return less(size() < rhs_size);
  }
  for (size_t i = rhs.size(); i > 0; --i) {
    if (digits_[i - 1 + shift] != rhs.digits_[i - 1]) {
      return less(digits_[i - 1 + shift] < rhs.digits_[i - 1]);
    }
  }
  for (size_t i = 0; i < shift; ++i) {
    if (digits_[i] != 0) {
      return std::strong_ordering::greater;
    }
  }
  return std::strong_ordering::equal;
}

bool big_integer::is_zero() const {
  return size() == 0;
}

big_integer& big_integer::abs_add_shifted(const big_integer& rhs, size_t shift) {
  if (rhs.is_zero()) {
    return *this;
  }
  digits_.resize(std::max(rhs.size() + shift, size()));
  bool carry = false;
  size_t i = 0;
  for (; i < rhs.size(); ++i) {
    digit& d = digits_[i + shift];
    digit rhs_d = rhs.digits_[i];
    digit sum = d + rhs_d + carry;
    carry = d > base - rhs_d || (d + rhs_d == base && carry);
    d = sum;
  }
  for (i += shift; i < size() && carry; ++i) {
    digits_[i]++;
    carry = digits_[i] == 0;
  }
  if (carry) {
    digits_.push_back(1);
  }
  check_invariant();
  return *this;
}

big_integer& big_integer::abs_sub_shifted(const big_integer& rhs, size_t shift) {
  digits_.resize(std::max(rhs.size() + shift, size()));
  bool borrow = false;
  bool is_smaller = std::strong_ordering::less == abs_compare_shifted(rhs, shift);
  if (is_smaller) {
    for (size_t i = 0; i < shift; ++i) {
      digits_[i] = base - digits_[i] + 1;
    }
    borrow = shift > 0;
  }
  size_t i = 0;
  for (; i < rhs.size(); ++i) {
    digit big_d = is_smaller ? rhs.digits_[i] : digits_[i + shift];
    digit small_d = is_smaller ? digits_[i + shift] : rhs.digits_[i];
    digit difference = big_d - small_d - borrow;
    borrow = (big_d == 0 && borrow) || big_d - borrow < small_d;
    digits_[i + shift] = difference;
  }
  for (i += shift; i < size() && borrow; ++i) {
    borrow = digits_[i] == 0;
    digits_[i]--;
  }
  is_negative_ = is_negative_ ^ is_smaller;
  strip_zeros();
  check_invariant();
  return *this;
}

big_integer& big_integer::add_shifted(const big_integer& rhs, size_t shift) {
  if (is_negative_ == rhs.is_negative_) {
    return abs_add_shifted(rhs, shift);
  }
  return abs_sub_shifted(rhs, shift);
}

big_integer& big_integer::sub_shifted(const big_integer& rhs, size_t shift) {
  return (this->negate().add_shifted(rhs, shift)).negate();
}

size_t big_integer::get_norm() const {
  digit d = digits_[size() - 1];
  size_t result = 0;
  while (d < base / 2 + 1) {
    d <<= 1;
    ++result;
  }
  return result;
}

std::pair<big_integer, big_integer> big_integer::divrem(const big_integer& rhs) const {
  if (size() < rhs.size()) {
    return {{}, *this};
  }
  int norm = static_cast<int>(rhs.get_norm());
  big_integer a = abs() << norm;
  big_integer b = rhs.abs() << norm;
  size_t n = b.size();
  size_t m = a.size() - n;
  big_integer quotient;
  bool is_big = a.abs_compare_shifted(b, m) != std::strong_ordering::less;
  quotient.digits_.resize(m + is_big);
  if (is_big) {
    a.abs_sub_shifted(b, m);
    quotient.digits_[m] = 1;
  }
  for (size_t i = m; i > 0; --i) {
    size_t j = i - 1;
    auto q = static_cast<digit>(
        ((static_cast<double_digit>(a.get_digit_value(n + j)) << digit_size) + a.get_digit_value(n + j - 1)) /
        b.digits_[n - 1]);
    quotient.digits_[j] = q;
    a.sub_shifted(b.mul_digit(q), j);
    while (a.is_negative_) {
      --quotient.digits_[j];
      a.add_shifted(b, j);
    }
  }
  return {quotient.negate_if(is_negative_ ^ rhs.is_negative_), 
          (a >> static_cast<int>(norm)).negate_if(is_negative_)};
}

big_integer& big_integer::to_abs() {
  is_negative_ = false;
  return *this;
}

big_integer big_integer::abs() const {
  return big_integer(*this).to_abs();
}

big_integer::big_integer(big_integer::digit d, bool is_negative) : digits_{d}, is_negative_(is_negative) {}

std::string to_string(const big_integer& a) {
  if (a.is_zero()) {
    return "0";
  }
  big_integer current(a);
  std::stack<big_integer::digit> decimal;
  while (!current.is_zero()) {
    auto [div, rem] = current.divrem(big_integer::base10);
    decimal.push(rem.is_zero() ? 0 : rem.digits_[0]);
    current = div;
  }
  std::string result = (a.is_negative_ ? "-" : "") + std::to_string(decimal.top());
  decimal.pop();
  while (!decimal.empty()) {
    std::string part = std::to_string(decimal.top());
    result += std::string(big_integer::exp10 - part.size(), '0') + part;
    decimal.pop();
  }
  return result;
}

big_integer& big_integer::to_zero() {
  digits_.clear();
  is_negative_ = false;
  return *this;
}

void big_integer::check_invariant() const {
  assert((size() == 0 && !is_negative_) || (digits_[size() - 1] != 0));
}

big_integer::digit big_integer::get_digit_value(size_t n) const {
  return n > size() ? 0 : digits_[n];
}

std::ostream& operator<<(std::ostream& out, const big_integer& a) {
  return out << to_string(a);
}
