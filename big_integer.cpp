#include "big_integer.h"

#include <algorithm>
#include <bit>
#include <cassert>
#include <cctype>
#include <cmath>
#include <compare>
#include <cstddef>
#include <functional>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>

#ifdef DEBUG
#define DEBUG_ONLY(expr) (expr)
#else
#define DEBUG_ONLY(_)
#endif

const size_t big_integer::DIGIT_SIZE = std::numeric_limits<big_integer::digit>::digits;
const big_integer::digit big_integer::MAX_DIGIT = std::numeric_limits<big_integer::digit>::max();
const big_integer::double_digit big_integer::BASE = static_cast<big_integer::double_digit>(big_integer::MAX_DIGIT) + 1;
const size_t big_integer::EXP10 = std::numeric_limits<digit>::digits10;
const big_integer big_integer::BASE10(static_cast<digit>(std::pow(10, big_integer::EXP10)), false);

static std::strong_ordering less(bool cond) {
  return cond ? std::strong_ordering::less : std::strong_ordering::greater;
}

static std::strong_ordering reverse(std::strong_ordering order) {
  return 0 <=> order;
}

big_integer::big_integer() : is_negative_{false} {}

big_integer::big_integer(const big_integer& other) = default;

big_integer::big_integer(const std::string& str) {
  if (str == "-0") {
    return;
  }
  is_negative_ = str.starts_with("-");
  if (str.empty() || std::any_of(str.begin() + is_negative_, str.end(), [](char c) { return c < 0 || !std::isdigit(c); })) {
    throw std::invalid_argument("Invalid string");
  }
  const size_t bound = EXP10 + is_negative_;
  size_t i = str.size();
  big_integer cur_base(1);
  while (true) {
    const bool in_bound = i >= bound;
    const size_t begin = in_bound ? i - EXP10 : is_negative_;
    const size_t end = in_bound ? EXP10 : i - is_negative_;
    abs_add_shifted(cur_base * big_integer(static_cast<digit>(std::stoul(str.substr(begin, end))), false));
    i = in_bound ? i - EXP10 : is_negative_;
    if (i == is_negative_) {
      break;
    }
    cur_base *= BASE10;
  }
  DEBUG_ONLY(check_invariant());
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
  size_t old_size = size();
  digits_.resize(size() + rhs.size());
  const auto begin = digits_.begin();
  const auto end = begin + static_cast<ptrdiff_t>(old_size);
  const auto shifted_begin = begin + static_cast<ptrdiff_t>(rhs.size());
  std::copy(begin, end, shifted_begin);
  std::fill(begin, shifted_begin, 0);
  for (size_t i = 0; i < old_size; ++i) {
    digit d = digits_[i + rhs.size()];
    digits_[i + rhs.size()] = 0;
    abs_add_shifted(rhs.mul_digit(d), i);
  }
  is_negative_ ^= rhs.is_negative_;
  strip_zeros();
  DEBUG_ONLY(check_invariant());
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
  const auto shift = static_cast<digit>(rhs);
  const size_t first = shift / DIGIT_SIZE;
  const size_t gain = first + (shift % DIGIT_SIZE != 0);
  const size_t remain = shift - first * DIGIT_SIZE;
  const size_t old_size = size();
  digits_.resize(old_size + gain);
  if (remain == 0) {
    std::shift_right(digits_.begin(), digits_.end(), static_cast<ptrdiff_t>(gain));
  } else {
    digit prev = 0;
    for (size_t i = old_size; i > 0; --i) {
      digits_[i - 1 + gain] = (prev << remain) | (digits_[i - 1] >> (DIGIT_SIZE - remain));
      prev = digits_[i - 1];
    }
    digits_[first] = digits_.front() << remain;
  }
  std::fill(digits_.begin(), digits_.begin() + static_cast<ptrdiff_t>(first), 0);
  strip_zeros();
  DEBUG_ONLY(check_invariant());
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
  const auto shift = static_cast<digit>(rhs);
  const size_t loss = shift / DIGIT_SIZE;
  if (loss >= size()) {
    return to_zero();
  }
  const size_t remain = shift - loss * DIGIT_SIZE;
  const bool is_neg = is_negative_;
  if (is_neg) {
    ++*this;
  }
  if (remain == 0) {
    std::shift_left(digits_.begin(), digits_.end(), static_cast<ptrdiff_t>(loss));
  } else {
    const size_t last = size() - loss - 1;
    for (size_t i = 0; i < last; ++i) {
      digits_[i] = (digits_[i + loss] >> remain) | (digits_[i + 1 + loss] << (DIGIT_SIZE - remain));
    }
    digits_[last] = digits_[last + loss] >> remain;
  }
  digits_.resize(size() - loss);
  strip_zeros();
  if (is_neg) {
    --*this;
  }
  DEBUG_ONLY(check_invariant());
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  big_integer result(*this);
  result.negate();
  return result;
}

big_integer big_integer::operator~() const {
  big_integer result(*this);
  --result.negate();
  return result;
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
  DEBUG_ONLY(check_invariant());
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
  big_integer result(a);
  result += b;
  return result;
}

big_integer operator-(const big_integer& a, const big_integer& b) {
  big_integer result(a);
  result -= b;
  return result;
}

big_integer operator*(const big_integer& a, const big_integer& b) {
  big_integer result(a);
  result *= b;
  return result;
}

big_integer operator/(const big_integer& a, const big_integer& b) {
  return a.divrem(b).first;
}

big_integer operator%(const big_integer& a, const big_integer& b) {
  return a.divrem(b).second;
}

big_integer operator&(const big_integer& a, const big_integer& b) {
  big_integer result(a);
  result &= b;
  return result;
}

big_integer operator|(const big_integer& a, const big_integer& b) {
  big_integer result(a);
  result |= b;
  return result;
}

big_integer operator^(const big_integer& a, const big_integer& b) {
  big_integer result(a);
  result ^= b;
  return result;
}

big_integer operator<<(const big_integer& a, int b) {
  big_integer result(a);
  result <<= b;
  return result;
}

big_integer operator>>(const big_integer& a, int b) {
  big_integer result(a);
  result >>= b;
  return result;
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

big_integer& big_integer::abs_add_digit(digit d) {
  if (d != 0) {
    digits_.resize(std::max(static_cast<size_t>(1), size()));
    bool carry = digits_.front() > MAX_DIGIT - d;
    digits_.front() += d;
    for (size_t i = 1; i < size() && carry; ++i) {
      digits_[i]++;
      carry = digits_[i] == 0;
    }
    if (carry) {
      digits_.push_back(1);
    }
  }
  return *this;
}

big_integer& big_integer::abs_sub_digit(digit d) {
  if (d != 0) {
    digits_.resize(std::max(static_cast<size_t>(1), size()));
    bool borrow = d > digits_.front();
    digits_.front() -= d;
    for (size_t i = 1; i < size() && borrow; ++i) {
      borrow = digits_[i] == 0;
      digits_[i]--;
    }
    if (size() == 1 && borrow) {
      digits_.front() = MAX_DIGIT - digits_.front() + 1;
      negate();
    } else {
      strip_zeros();
    }
  }
  return *this;
}

big_integer& big_integer::add_digit(digit d) {
  if (is_negative_) {
    abs_sub_digit(d);
  } else {
    abs_add_digit(d);
  }
  return *this;
}

big_integer& big_integer::sub_digit(digit d) {
  if (is_negative_) {
    abs_add_digit(d);
  } else {
    abs_sub_digit(d);
  }
  return *this;
}

big_integer& big_integer::negate() {
  is_negative_ = !is_zero() && !is_negative_;
  DEBUG_ONLY(check_invariant());
  return *this;
}

big_integer& big_integer::negate_if(bool cond) {
  is_negative_ = cond && !is_zero() && !is_negative_;
  DEBUG_ONLY(check_invariant());
  return *this;
}

big_integer& big_integer::bitwise(const big_integer& rhs,
                                  const std::function<big_integer::digit(big_integer::digit, big_integer::digit)>& f) {
  if (size() < rhs.size()) {
    digits_.resize(rhs.size());
  }
  bool borrow = is_negative_, rhs_borrow = rhs.is_negative_;
  const bool is_neg = f(is_negative_, rhs.is_negative_);
  size_t i = 0;
  for (; i < rhs.size(); ++i) {
    const digit d = digits_[i] - borrow;
    const digit rhs_d = rhs.digits_[i] - rhs_borrow;
    borrow = digits_[i] == 0 && borrow;
    rhs_borrow = rhs.digits_[i] == 0 && rhs_borrow;
    digits_[i] = f(is_negative_ ? ~d : d, rhs.is_negative_ ? ~rhs_d : rhs_d) ^ (is_neg ? MAX_DIGIT : 0);
  }
  for (; i < size(); ++i) {
    digit d = digits_[i] - borrow;
    digits_[i] = f(is_negative_ ? ~d : d, (rhs.is_negative_ ? MAX_DIGIT : 0)) ^ (is_neg ? MAX_DIGIT : 0);
  }
  is_negative_ = is_neg;
  strip_zeros();
  if (is_neg) {
    --*this;
  }
  DEBUG_ONLY(check_invariant());
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
    const double_digit product = static_cast<double_digit>(digits_[i]) * d;
    const digit lo = product;
    digit& rd = result.digits_[i];
    const digit sum = rd + product + carry;
    carry = rd > MAX_DIGIT - lo || (rd + lo == MAX_DIGIT && carry);
    rd = sum;
    result.digits_[i + 1] += product >> DIGIT_SIZE;
  }
  const double_digit product = static_cast<double_digit>(digits_.back()) * d;
  const digit lo = product;
  digit& rd = result.digits_.back();
  const digit sum = rd + product + carry;
  carry = rd > MAX_DIGIT - lo || (rd + lo == MAX_DIGIT && carry);
  rd = sum;
  digit hi = product >> DIGIT_SIZE;
  if (hi != 0 || carry) {
    result.digits_.push_back(hi + carry);
    if (hi == MAX_DIGIT && carry) {
      result.digits_.push_back(1);
    }
  }
  DEBUG_ONLY(result.check_invariant());
  return result;
}

size_t big_integer::size() const {
  return digits_.size();
}

void big_integer::strip_zeros() {
  while (!digits_.empty() && digits_.back() == 0) {
    digits_.pop_back();
  }
  is_negative_ &= !is_zero();
  DEBUG_ONLY(check_invariant());
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
    const digit rhs_d = rhs.digits_[i];
    const digit sum = d + rhs_d + carry;
    carry = d > MAX_DIGIT - rhs_d || (d + rhs_d == MAX_DIGIT && carry);
    d = sum;
  }
  for (i += shift; i < size() && carry; ++i) {
    digits_[i]++;
    carry = digits_[i] == 0;
  }
  if (carry) {
    digits_.push_back(1);
  }
  return *this;
}

big_integer& big_integer::abs_sub_shifted(const big_integer& rhs, size_t shift) {
  digits_.resize(std::max(rhs.size() + shift, size()));
  bool borrow = false;
  const bool is_smaller = std::strong_ordering::less == abs_compare_shifted(rhs, shift);
  if (is_smaller) {
    for (size_t i = 0; i < shift; ++i) {
      digits_[i] = MAX_DIGIT - digits_[i] + 1;
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
  DEBUG_ONLY(check_invariant());
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
  return std::max(static_cast<int>(DIGIT_SIZE) - static_cast<int>(std::bit_width(digits_.back())),
                  static_cast<int>(0));
}

std::pair<big_integer, big_integer> big_integer::divrem(const big_integer& rhs) const {
  if (size() < rhs.size()) {
    return {{}, *this};
  }
  const int norm = static_cast<int>(rhs.get_norm());
  big_integer a = abs() << norm;
  const big_integer b = rhs.abs() << norm;
  const size_t n = b.size();
  const size_t m = a.size() - n;
  big_integer quotient;
  const bool is_big = a.abs_compare_shifted(b, m) != std::strong_ordering::less;
  quotient.digits_.resize(m + is_big);
  if (is_big) {
    a.abs_sub_shifted(b, m);
    quotient.digits_[m] = 1;
  }
  for (size_t i = m; i > 0; --i) {
    const size_t j = i - 1;
    const auto q = static_cast<digit>(
        ((static_cast<double_digit>(a.get_digit_value(n + j)) << DIGIT_SIZE) + a.get_digit_value(n + j - 1)) /
        b.digits_[n - 1]);
    quotient.digits_[j] = q;
    a.sub_shifted(b.mul_digit(q), j);
    while (a.is_negative_) {
      --quotient.digits_[j];
      a.add_shifted(b, j);
    }
  }
  return {quotient.negate_if(is_negative_ ^ rhs.is_negative_), (a >> static_cast<int>(norm)).negate_if(is_negative_)};
}

big_integer& big_integer::to_abs() {
  is_negative_ = false;
  return *this;
}

big_integer big_integer::abs() const {
  big_integer result(*this);
  result.to_abs();
  return result;
}

big_integer::big_integer(big_integer::digit d, bool is_negative) : digits_{d}, is_negative_(is_negative) {}

std::string to_string(const big_integer& a) {
  if (a.is_zero()) {
    return "0";
  }
  big_integer current(a);
  std::vector<big_integer::digit> decimal;
  while (!current.is_zero()) {
    const auto [div, rem] = current.divrem(big_integer::BASE10);
    decimal.push_back(rem.is_zero() ? 0 : rem.digits_.front());
    current = div;
  }
  std::string result = (a.is_negative_ ? "-" : "") + std::to_string(decimal.back());
  decimal.pop_back();
  while (!decimal.empty()) {
    const std::string part = std::to_string(decimal.back());
    result += std::string(big_integer::EXP10 - part.size(), '0') + part;
    decimal.pop_back();
  }
  return result;
}

big_integer& big_integer::to_zero() {
  digits_.clear();
  is_negative_ = false;
  return *this;
}

void big_integer::check_invariant() const {
  assert((size() == 0 && !is_negative_) || (digits_.back() != 0));
}

big_integer::digit big_integer::get_digit_value(size_t n) const {
  return n >= size() ? 0 : digits_[n];
}

std::ostream& operator<<(std::ostream& out, const big_integer& a) {
  return out << to_string(a);
}
