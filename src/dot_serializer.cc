#include "lancet2/dot_serializer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>

#include "absl/strings/str_format.h"
#include "lancet2/core_enums.h"
#include "lancet2/utils.h"

namespace lancet2 {
void Graph::DotSerializer::WriteComponent(usize comp_id, const std::string& suffix) const {
  const auto outDir = graphPtr->params->outGraphsDir.empty() ? std::filesystem::current_path()
                                                             : std::filesystem::path(graphPtr->params->outGraphsDir);
  const auto windowID = graphPtr->window->ToRegionString();
  const auto fName = absl::StrFormat("%s/%s_c%d_%s.dot", outDir, windowID, comp_id, suffix);
  std::ofstream outStream(fName, std::ios::trunc);

  DumpHeader(outStream);
  DumpComponent(comp_id, outStream);
  outStream << "}\n";
  outStream.close();
}

void Graph::DotSerializer::WriteComponent(usize comp_id, absl::Span<const PathNodeIds> flow_paths) const {
  //  https://www.wolframalpha.com/input/?i=Golden+ratio
  static constexpr double goldenRatioConjugate = 0.618033988749894848204586834365638117720309179805762862135;
  usize pathNum = 0;
  double currentHue = 0.0;

  const auto windowID = graphPtr->window->ToRegionString();
  const auto outDir = graphPtr->params->outGraphsDir.empty() ? std::filesystem::current_path()
                                                             : std::filesystem::path(graphPtr->params->outGraphsDir);
  const auto fName = absl::StrFormat("%s/%s_c%d_path_flow.dot", outDir, windowID, comp_id);
  std::ofstream outStream(fName, std::ios::out | std::ios::trunc);
  DumpHeader(outStream);
  DumpComponent(comp_id, outStream);

  std::for_each(flow_paths.cbegin(), flow_paths.cend(), [&](const PathNodeIds& path_flow) {
    pathNum++;
    currentHue += goldenRatioConjugate;
    currentHue = std::fmod(currentHue, 1.0);
    DumpPathFlow(path_flow, currentHue, outStream);
  });

  outStream << "}\n";
  outStream.close();
}

void Graph::DotSerializer::DumpHeader(std::ostream& out_stream) {
  static constexpr auto header = R"raw(strict digraph G {
graph [layout=neato,bgcolor=black,size="120,180",ratio=compress,rankdir=LR,overlap=vpsc,overlap_shrink=true,start=self];
node [style=filled,fontsize=2,width=2,height=2,fixedsize=false];
edge [color=gray,fontsize=8,fontcolor=floralwhite,len=3,fixedsize=false,headclip=true,tailclip=true];
)raw";

  out_stream << header;
}

void Graph::DotSerializer::DumpComponent(usize comp_id, std::ostream& out_stream) const {
  out_stream << absl::StreamFormat("subgraph component_%d", comp_id) << " {\n";

  std::for_each(graphPtr->nodesMap.cbegin(), graphPtr->nodesMap.cend(), [&](Graph::NodeContainer::const_reference p) {
    if (comp_id != 0 && p.second->ComponentID != comp_id) return;

    const auto isSrc = p.second->IsSource();
    const auto isSnk = p.second->IsSink();

    const auto forwardSequence = p.second->GetFwdSeq();
    const auto nodeName = isSrc ? "source" : isSnk ? "sink" : absl::StrFormat("nodeID_%d", p.first);
    const auto* nodeShape = isSrc || isSnk ? "diamond" : "circle";

    out_stream << absl::StreamFormat(R"raw(%d [shape=%s fillcolor=%s label="F:%s\nR:%s\n %s:%s\ncoverage=%d"]
)raw",
                                     p.first, nodeShape, p.second->FillColor(), forwardSequence,
                                     OppositeStrandSequence(forwardSequence), nodeName,
                                     ToString(p.second->GetOrientation()), p.second->TotalSampleCount());

    for (const Edge& e : *p.second) {
      out_stream << absl::StreamFormat("%d -> %d [taillabel=%s headlabel=%s]\n", p.first, e.GetDstID(),
                                       ToString(e.GetSrcDir()), ToString(e.GetDstDir()));
    }
  });

  out_stream << "}\n";
}

void Graph::DotSerializer::DumpPathFlow(const PathNodeIds& path_flow, double hue, std::ostream& out_stream) {
  std::for_each(path_flow.cbegin(), path_flow.cend(), [&out_stream, &hue](const EdgeNodeIds& ep) {
    out_stream << absl::StreamFormat("%d -> %d [color=\"%f 0.85 0.85\"]", ep.srcId, ep.dstId, hue);
  });
}

auto Graph::DotSerializer::OppositeStrandSequence(absl::string_view seq) -> std::string {
  std::string result(seq.cbegin(), seq.cend());
  std::transform(result.begin(), result.end(), result.begin(), [](const char& b) { return utils::RevComp(b); });
  return result;
}
}  // namespace lancet2
