#include "big_integer.h"

#include <cstddef>
#include <cstring>
#include <functional>
#include <limits>
#include <ostream>
#include <stdexcept>


static bool are_zeros(std::vector<big_integer::digit>::const_iterator begin,
                      std::vector<big_integer::digit>::const_iterator end) {
  return std::all_of(begin, end, [](big_integer::digit d) { return d == 0; });
}

static std::strong_ordering less(bool cond) { 
  return cond ? std::strong_ordering::less : std::strong_ordering::greater;
}

big_integer::big_integer() : is_negative_{false} {}

big_integer::big_integer(const big_integer& other) = default;

big_integer::big_integer(int a) 
  : digits_{a == std::numeric_limits<int>::min() ? 
      static_cast<digit>(std::numeric_limits<int>::max()) + 1 :
      static_cast<digit>(std::abs(a))}
  , is_negative_{a < 0} {}

big_integer::big_integer(const std::string& str) {}

big_integer::~big_integer() = default;

big_integer& big_integer::operator=(const big_integer& other) = default;

big_integer& big_integer::operator+=(const big_integer& rhs) {
  digits_.resize(std::max(rhs.size(), size()));
  if (is_negative_ == rhs.is_negative_) {
    bool carry = false;
    size_t i = 0;
    for (; i < rhs.size(); ++i) {
      digit sum = digits_[i] + rhs.digits_[i] + carry;
      carry = digits_[i] > big_integer::base - rhs.digits_[i] || (digits_[i] + rhs.digits_[i] == big_integer::base && carry);
      digits_[i] = sum;
    }
    for (; i < size() && carry; ++i) {
      digits_[i]++;
      carry = digits_[i] == 0;
    }
    if (carry) {
      digits_.push_back(1);
    }
    
  } else {
    bool all_zero = true;
    bool is_smaller = (negate() < rhs) ^ rhs.is_negative_;
    negate();
    const big_integer &smaller = is_smaller ? *this : rhs,
                      &bigger  = is_smaller ? rhs : *this;
    bool borrow = false;
    size_t i = 0;
    for (; i < smaller.size(); ++i) {
      digit difference = bigger.digits_[i] - smaller.digits_[i] - borrow;
      borrow = bigger.digits_[i] == 0 || bigger.digits_[i] - borrow < smaller.digits_[i];
      digits_[i] = difference;
      all_zero &= digits_[i] == 0;
    }
    for (; i < size() && borrow; ++i) {
      borrow = digits_[i] == 0;
      digits_[i]--;
      all_zero &= digits_[i] == 0;
    }
    is_negative_ = !all_zero && (is_negative_ ^ is_smaller);
  }
  return *this;
}

big_integer& big_integer::operator-=(const big_integer& rhs) {
  return (this->negate() += rhs).negate();
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
  digit shift = static_cast<digit>(rhs);
  size_t start = shift / big_integer::digit_size;
  size_t gain = start + (shift % big_integer::digit_size != 0);
  size_t remain = shift - start * big_integer::digit_size;
  size_t old_size = size();
  digits_.resize(old_size + gain);
  digit prev = 0;
  for (size_t i = old_size; i > 0; --i) {
    digits_[i - 1 + gain] = (prev << remain) | (digits_[i - 1] >> (big_integer::digit_size - remain));
    prev = digits_[i - 1];
  }
  digits_[start] = digits_[0] << remain;
  std::fill(digits_.begin(), digits_.begin() + static_cast<ptrdiff_t>(start), 0);
  if (is_negative_) {
    --*this;
  }
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
  digit shift = static_cast<digit>(rhs);
  size_t loss = shift / big_integer::digit_size;
  size_t remain = shift - loss * big_integer::digit_size;
  digit prev = 0;
  size_t last = size() - loss - 1;
  for (size_t i = 0; i < last; ++i) {
    digits_[i] = (digits_[i + loss] >> remain) | (digits_[i + 1 + loss] << (big_integer::digit_size - remain));
  }
  digits_[last] = digits_[last + loss] >> remain;
  std::fill(digits_.end() - static_cast<ptrdiff_t>(loss), digits_.end(), 0);
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
    bool all_zero = true;
    bool borrow = true;
    for (size_t i = 0; i < size() && borrow; ++i) {
      borrow = digits_[i] == 0;
      digits_[i]--;
      all_zero = digits_[i] == 0;
    }
    is_negative_ &= !all_zero;
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
  result.is_negative_ = a.is_negative_;
  for (size_t i = 0; i < b.size(); ++i) {
    result += a.mul_digit(b.digits_[i]) << static_cast<int>(big_integer::digit_size * i);
  }
  result.is_negative_ = a.is_negative_ ^ b.is_negative_;
  return result;
}

big_integer operator/(const big_integer& a, const big_integer& b) {

  return a;
}

big_integer operator%(const big_integer& a, const big_integer& b) {
  return a;
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
  return big_integer(a) <<= b;
}

std::strong_ordering operator<=>(const big_integer& a, const big_integer& b) {
  if (a.is_negative_ != b.is_negative_) {
    return less(a.is_negative_);
  }
  size_t common_size = std::min(a.size(), b.size());
  if (!are_zeros(a.digits_.begin() + static_cast<ptrdiff_t>(common_size), a.digits_.end())) {
    return less(a.is_negative_);
  } else if (!are_zeros(b.digits_.begin() + static_cast<ptrdiff_t>(common_size), b.digits_.end())) {
    return less(!a.is_negative_);
  }
  for (size_t i = common_size; i > 0; --i) {
    if (a.digits_[i - 1] != b.digits_[i - 1]) {
      return less((a.digits_[i - 1] < b.digits_[i - 1]) ^ a.is_negative_);
    }
  }
  return std::strong_ordering::equal;
}

bool operator==(const big_integer& a, const big_integer& b) = default;
bool operator!=(const big_integer& a, const big_integer& b) = default;
bool operator<(const big_integer& a, const big_integer& b) = default;
bool operator>(const big_integer& a, const big_integer& b) = default;
bool operator<=(const big_integer& a, const big_integer& b) = default;
bool operator>=(const big_integer& a, const big_integer& b) = default;

big_integer& big_integer::negate() {
  is_negative_ = !is_negative_;
  return *this;
}

big_integer& big_integer::bitwise(const big_integer& rhs, std::function<big_integer::digit(big_integer::digit, big_integer::digit)> f) {
  size_t common_size = std::min(size(), rhs.size());
  bool borrow = true, rhs_borrow = true;
  bool is_neg = f(is_negative_, rhs.is_negative_);
  size_t i = 0;
  for (; i < common_size; ++i) {
    borrow = digits_[i] == 0 && borrow;
    rhs_borrow = rhs.digits_[i] == 0 && rhs_borrow;
    digit d = digits_[i] - (is_negative_ && borrow ? 1 : 0);
    digit rhs_d = rhs.digits_[i] - (rhs.is_negative_ && rhs_borrow ? 1 : 0);
    digits_[i] = f((d ^ (is_negative_ ? big_integer::base : 0)), (rhs_d ^ (rhs.is_negative_ ? big_integer::base : 0))) ^ (is_neg ? big_integer::base : 0);
    digits_[i] ^= (is_neg ? big_integer::base : 0);
  }
  for (; i < digits_.size(); ++i) {
    digit d = digits_[i] - (is_negative_ && borrow ? 1 : 0);
    digits_[i] = f((d ^ (is_negative_ ? big_integer::base : 0)), (rhs.is_negative_ ? base : 0)) ^ (is_neg ? big_integer::base : 0);
  }

  is_negative_ = is_neg;
  if (is_neg) {
    ++*this;
  }
  return *this;
}

big_integer big_integer::mul_digit(digit d) const {
  big_integer result;
  result.digits_.resize(size());
  for (size_t i = 0; i < size() - 1; ++i) {
    double_digit product = static_cast<double_digit>(digits_[i]) * static_cast<double_digit>(d);
    result.digits_[i] += product;
    result.digits_[i + 1] += product / digit_size;
  }
  double_digit product = static_cast<double_digit>(digits_[size() - 1]) * d;
  result.digits_[size() - 1] += product;
  digit hi = product / digit_size;
  if (hi != 0) {
    result.digits_.push_back(hi);
  }
  return result;
}

size_t big_integer::size() const {
  return digits_.size();
}

std::string to_string(const big_integer& a) {
  return "";
}

std::ostream& operator<<(std::ostream& out, const big_integer& a) {
  return out << to_string(a);
}


