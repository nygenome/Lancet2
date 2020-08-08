#include "lancet/vcf_writer.h"

#include <stdexcept>

#include "absl/strings/str_format.h"
#include "htslib/bgzf.h"
#include "htslib/tbx.h"

namespace lancet {
class VcfWriter::Impl {
 public:
  explicit Impl(const std::filesystem::path& out_path) : vcfPath(out_path) { OpenBgzfHandle(out_path.c_str()); }

  auto Write(std::string_view record) -> absl::Status {
    const auto numWritten = bgzf_write(fp, record.data(), record.length());
    return numWritten == record.length() ? absl::OkStatus() : absl::InternalError("could not write to BGZF handle");
  }

  void Flush() {
    if (fp == nullptr || isClosed) return;
    if (bgzf_flush(fp) != 0) throw std::runtime_error("could not flush variants to vcf");
    return;
  }

  void Close() {
    if (fp == nullptr || isClosed) return;
    bgzf_close(fp);
    tbx_index_build(vcfPath.c_str(), 0, &tbx_conf_vcf);
  }

 private:
  bool isClosed = false;
  std::filesystem::path vcfPath;
  BGZF* fp = nullptr;

  void OpenBgzfHandle(const char* path) {
    fp = bgzf_open(path, "w");
    if (fp == nullptr) {
      const auto errMsg = absl::StrFormat("could not open BGZF handle for %s", path);
      throw std::runtime_error(errMsg);
    }
  }
};

VcfWriter::VcfWriter(const std::filesystem::path& out_path) : pimpl(std::make_unique<Impl>(out_path)) {}
VcfWriter::~VcfWriter() { Close(); }
auto VcfWriter::Write(std::string_view record) -> absl::Status { return pimpl->Write(record); }
void VcfWriter::Flush() { return pimpl->Flush(); }
void VcfWriter::Close() { return pimpl->Close(); }
}  // namespace lancet
