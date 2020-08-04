/* Copyright 2017 The TensorFlow Authors. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#pragma once

#include "absl/status/status.h"

namespace lancet::internal_statusor {

class Helper {
 public:
  // Move type-agnostic error handling to the .cc.
  static void HandleInvalidStatusCtorArg(absl::Status*);
  [[noreturn]] static void Crash(const absl::Status& status);
};

// Construct an instance of T in `p` through placement new, passing Args... to
// the constructor.
// This abstraction is here mostly for the gcc performance fix.
template <typename T, typename... Args>
void PlacementNew(void* p, Args&&... args) {
#if defined(__GNUC__) && !defined(__clang__)
  // Teach gcc that 'p' cannot be null, fixing code size issues.
  if (p == nullptr) __builtin_unreachable();
#endif
  new (p) T(std::forward<Args>(args)...);
}

// Helper base class to hold the data and all operations.
// We move all this to a base class to allow mixing with the appropriate
// TraitsBase specialization.
template <typename T>
class StatusOrData {
  template <typename U>
  friend class StatusOrData;

 public:
  StatusOrData() = delete;

  StatusOrData(const StatusOrData& other) {
    if (other.ok()) {
      MakeValue(other.data_);
      MakeStatus();
    } else {
      MakeStatus(other.status_);
    }
  }

  StatusOrData(StatusOrData&& other) noexcept {
    if (other.ok()) {
      MakeValue(std::move(other.data_));
      MakeStatus();
    } else {
      MakeStatus(other.status_);
    }
  }

  template <typename U>
  StatusOrData(const StatusOrData<U>& other) {  // NOLINT
    if (other.ok()) {
      MakeValue(other.data_);
      MakeStatus();
    } else {
      MakeStatus(other.status_);
    }
  }

  template <typename U>
  StatusOrData(StatusOrData<U>&& other) {  // NOLINT
    if (other.ok()) {
      MakeValue(std::move(other.data_));
      MakeStatus();
    } else {
      MakeStatus(other.status_);
    }
  }

  explicit StatusOrData(const T& value) : data_(value) { MakeStatus(); }
  explicit StatusOrData(T&& value) : data_(std::move(value)) { MakeStatus(); }

  explicit StatusOrData(const absl::Status& status) : status_(status) { EnsureNotOk(); }  // NOLINT
  explicit StatusOrData(absl::Status&& status) : status_(std::move(status)) { EnsureNotOk(); }

  auto operator=(const StatusOrData& other) -> StatusOrData& {
    if (this == &other) return *this;
    if (other.ok()) {
      Assign(other.data_);
    } else {
      Assign(other.status_);
    }
    return *this;
  }

  auto operator=(StatusOrData&& other) noexcept -> StatusOrData& {
    if (this == &other) return *this;
    if (other.ok()) {
      Assign(std::move(other.data_));
    } else {
      Assign(std::move(other.status_));
    }
    return *this;
  }

  ~StatusOrData() {
    if (ok()) {
      status_.~Status();  // NOLINT
      data_.~T();         // NOLINT
    } else {
      status_.~Status();  // NOLINT
    }
  }

  void Assign(const T& value) {
    if (ok()) {
      data_.~T();  // NOLINT
      MakeValue(value);
    } else {
      MakeValue(value);
      status_ = absl::OkStatus();  // NOLINT
    }
  }

  void Assign(T&& value) {
    if (ok()) {
      data_.~T();  // NOLINT
      MakeValue(std::move(value));
    } else {
      MakeValue(std::move(value));
      status_ = absl::OkStatus();  // NOLINT
    }
  }

  void Assign(const absl::Status& status) {
    Clear();
    status_ = status;  // NOLINT
    EnsureNotOk();
  }

  void Assign(absl::Status&& status) {
    Clear();
    // Note that we copy instead of moving the status here so that
    // status.~StatusOrData() can call ok() without invoking UB.
    status_ = status;  // NOLINT
    EnsureNotOk();
  }

  [[nodiscard]] auto ok() const -> bool {
    return status_.ok();  // NOLINT
  }

 protected:
  // status_ will always be active after the constructor.
  // We make it a union to be able to initialize exactly how we need without
  // waste.
  // Eg. in the copy constructor we use the default constructor of Status in
  // the ok() path to avoid an extra Ref call.
  union {  // NOLINT
    absl::Status status_;
  };

  // data_ is active iff status_.ok()==true
  struct Dummy {};
  union {  // NOLINT
    // When T is const, we need some non-const object we can cast to void* for
    // the placement new. dummy_ is that object.
    Dummy dummy_;
    T data_;
  };

  void Clear() {
    if (ok()) data_.~T();  // NOLINT
  }

  void EnsureOk() const {
    if (!ok()) Helper::Crash(status_);  // NOLINT
  }

  void EnsureNotOk() {
    if (ok()) Helper::HandleInvalidStatusCtorArg(&status_);  // NOLINT
  }

  // Construct the value (ie. data_) through placement new with the passed
  // argument.
  template <typename Arg>
  void MakeValue(Arg&& arg) {
    internal_statusor::PlacementNew<T>(&dummy_, std::forward<Arg>(arg));  // NOLINT
  }

  // Construct the status (ie. status_) through placement new with the passed
  // argument.
  template <typename... Args>
  void MakeStatus(Args&&... args) {
    internal_statusor::PlacementNew<absl::Status>(&status_, std::forward<Args>(args)...);  // NOLINT
  }
};

// Helper base class to allow implicitly deleted constructors and assignment
// operations in StatusOr.
// TraitsBase will explicitly delete what it can't support and StatusOr will
// inherit that behavior implicitly.
template <bool Copy, bool Move>
struct TraitsBase {
  TraitsBase() = default;
  ~TraitsBase() = default;
  TraitsBase(const TraitsBase&) = default;
  TraitsBase(TraitsBase&&) noexcept = default;
  auto operator=(const TraitsBase&) -> TraitsBase& = default;
  auto operator=(TraitsBase&&) noexcept -> TraitsBase& = default;
};

template <>
struct TraitsBase<false, true> {
  TraitsBase() = default;
  ~TraitsBase() = default;
  TraitsBase(const TraitsBase&) = delete;
  TraitsBase(TraitsBase&&) = default;
  auto operator=(const TraitsBase&) -> TraitsBase& = delete;
  auto operator=(TraitsBase&&) noexcept -> TraitsBase& = default;
};

template <>
struct TraitsBase<false, false> {
  TraitsBase() = default;
  ~TraitsBase() = default;
  TraitsBase(const TraitsBase&) = delete;
  TraitsBase(TraitsBase&&) = delete;
  auto operator=(const TraitsBase&) -> TraitsBase& = delete;
  auto operator=(TraitsBase&&) noexcept -> TraitsBase& = delete;
};

}  // namespace lancet::internal_statusor
