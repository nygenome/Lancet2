#ifndef SRC_LANCET_CORE_WINDOW_H_
#define SRC_LANCET_CORE_WINDOW_H_

#include "lancet/base/types.h"
#include "lancet/hts/reference.h"

#include "spdlog/fmt/bundled/core.h"

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

namespace lancet::core {

class Window {
 public:
  using Chrom = hts::Reference::Chrom;
  using RefPath = std::filesystem::path;
  using RegSpec = hts::Reference::ParseRegionResult;
  using RegionPtr = std::shared_ptr<hts::Reference::Region const>;

  Window() = default;
  Window(RegSpec reg_spec, Chrom chrom, RefPath ref_path)
      : mSpec(std::move(reg_spec)), mChrom(std::move(chrom)), mRefPath(std::move(ref_path)) {}

  void SetGenomeIndex(usize const window_index) { mGenIdx = window_index; }

  [[nodiscard]] auto GenomeIndex() const -> usize { return mGenIdx; }
  [[nodiscard]] auto ChromIndex() const -> usize { return mChrom.Index(); }
  [[nodiscard]] auto ChromName() const -> std::string { return mChrom.Name(); }
  [[nodiscard]] auto StartPos1() const -> u64 { return mSpec.mRegionSpan[0].value_or(1); }
  [[nodiscard]] auto EndPos1() const -> u64 {
    return mSpec.mRegionSpan[1].value_or(mChrom.Length());
  }
  [[nodiscard]] auto Length() const -> usize {
    return mSpec.Length() != 0 ? mSpec.Length() : mChrom.Length();
  }

  [[nodiscard]] auto ToSamtoolsRegion() const -> std::string {
    auto const name_has_colon = mSpec.mChromName.find(':') != std::string::npos;
    return name_has_colon ? fmt::format("{{{}}}:{}-{}", ChromName(), StartPos1(), EndPos1())
                          : fmt::format("{}:{}-{}", ChromName(), StartPos1(), EndPos1());
  }

  [[nodiscard]] auto SeqView() const -> std::string_view {
    EnsureRegionBuilt();
    return mRegPtr->SeqView();
  }

  [[nodiscard]] auto SeqData() const -> char const* {
    EnsureRegionBuilt();
    return mRegPtr->SeqData();
  }

  [[nodiscard]] auto AsRegionPtr() const -> RegionPtr {
    EnsureRegionBuilt();
    return mRegPtr;
  }

 private:
  usize mGenIdx = 0;
  Chrom mChrom;
  RegSpec mSpec;
  RefPath mRefPath;
  mutable RegionPtr mRegPtr = nullptr;

  void EnsureRegionBuilt() const {
    if (mRegPtr != nullptr || mRefPath.empty() || mSpec.mChromName.empty())
      return;

    hts::Reference const reference(mRefPath);
    mRegPtr = std::make_shared<hts::Reference::Region const>(reference.MakeRegion(mSpec));
  }
};

using WindowPtr = std::shared_ptr<Window>;

}  // namespace lancet::core

#endif  // SRC_LANCET_CORE_WINDOW_H_
