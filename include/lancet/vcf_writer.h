#pragma once

#include <filesystem>
#include <memory>
#include <string_view>

#include "absl/status/status.h"

namespace lancet {
class VcfWriter {
 public:
  explicit VcfWriter(const std::filesystem::path& out_path);
  VcfWriter() = delete;
  ~VcfWriter();

  VcfWriter(VcfWriter&&) noexcept = default;
  auto operator=(VcfWriter&&) noexcept -> VcfWriter& = default;
  VcfWriter(const VcfWriter&) = delete;
  auto operator=(const VcfWriter&) -> VcfWriter& = delete;

  auto Write(std::string_view rec) -> absl::Status;

  void Flush();
  void Close();

 private:
  class Impl;
  std::unique_ptr<Impl> pimpl;
};
}  // namespace lancet
