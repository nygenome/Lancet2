#include "lancet2/gfa_serializer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>

#include "absl/strings/str_split.h"
#include "absl/strings/str_format.h"
#include "lancet2/core_enums.h"
#include "lancet2/utils.h"

namespace lancet2 {

void Graph::GfaSerializer::WriteComponent(usize comp_id, const std::string& suffix) const { 
  const auto outDir = graphPtr->params->outGraphsDir.empty() ? std::filesystem::current_path() :
                                                               std::filesystem::path(graphPtr->params->outGraphsDir);
  const auto windowID = graphPtr->window->ToRegionString();
  const auto fName = absl::StrFormat("%s/%s_c%d_%s.gfa", outDir, windowID, comp_id, suffix);
  std::ofstream outStream(fName, std::ios::trunc);

  DumpHeader(outStream);
  DumpComponent(comp_id, outStream);
  outStream.close();
}


void Graph::GfaSerializer::WriteComponent(usize comp_id, absl::Span<const PathNodeIds> flow_paths, std::string& windowId) const {
  // https:://www.wolframalpha.com/input/?i=Golden+ratio
  static constexpr double goldenRatioConjugate = 0.618033988749894848204586834365638117720309179805762862135;
  usize pathNum = 0;
  double currentHue = 0.0;

  const auto windowID = graphPtr->window->ToRegionString();
  const auto outDir = graphPtr->params->outGraphsDir.empty() ? std::filesystem::current_path()
                                                             : std::filesystem::path(graphPtr->params->outGraphsDir);
  const auto fName = absl::StrFormat("%s/%s_c%d_path_flow.gfa", outDir, windowID, comp_id);
  std::ofstream outStream(fName, std::ios::out | std::ios::trunc);
  DumpHeader(outStream);
  DumpComponent(comp_id, outStream);
  
  std::for_each(flow_paths.cbegin(), flow_paths.cend(), [&](const PathNodeIds& path_flow) {
    pathNum++;
    
    DumpPathFlow(path_flow, pathNum, windowId, outStream);
  }); 
  outStream.close();
}

void Graph::GfaSerializer::DumpHeader(std::ostream& out_stream) {
  static constexpr auto header = "H\tVN:Z:1.0\n";
  out_stream << header;
}

void Graph::GfaSerializer::DumpComponent(usize comp_id, std::ostream& out_stream) const {
  std::ostringstream linkStream;
  std::for_each(graphPtr->nodesMap.cbegin(), graphPtr->nodesMap.cend(), [&](Graph::NodeContainer::const_reference p) {
    if (comp_id != 0 && p.second->ComponentID != comp_id) return;

    const auto isSrc = p.second->IsSource();
    const auto isSnk = p.second->IsSink();

    const auto forwardSequence = p.second->GetFwdSeq();
    const auto nodeName = isSrc ? "source" : isSnk ? "sink" : absl::StrFormat("nodeID_%d", p.first);
    const auto* nodeShape = isSrc || isSnk ? "diamond" : "circle";

    auto node_label = p.second->GetLabel();
    const auto hasRef = node_label.HasLabel(KmerLabel::REFERENCE);
    const auto hasTmr = node_label.HasLabel(KmerLabel::TUMOR);
    const auto hasNml = node_label.HasLabel(KmerLabel::NORMAL);
    auto final_label = "";
    if (hasRef || (hasTmr && hasNml)) final_label = "REF/SHARED";
    if (hasTmr && !hasRef && !hasNml) final_label = "TUMOR";
    if (hasNml && !hasRef && !hasTmr) final_label = "NORMAL";
    out_stream << absl::StreamFormat("S\t%d\t%s\tLN:i:%d\ttc:i:%d\tnc:i:%d\tor:i:%d\tla:Z:%s\n",p.first,forwardSequence,p.second->GetLength(),p.second->SampleCount(SampleLabel::TUMOR),p.second->SampleCount(SampleLabel::NORMAL),p.second->GetOrientation(),final_label);

    for (const Edge& e : *p.second) {
      linkStream << absl::StreamFormat("L\t%d\t+\t%d\t+\t%dM\n",p.first,e.GetDstID(),(p.second->GetLength() - 1));
      //out_stream << absl::StreamFormat("L\t%d\t+\t%d\t+\t%dM\n",p.first, e.GetDstID(),(p.second->GetLength() - 1));
                                       //ToString(e.GetSrcDir()), ToString(e.GetDstDir()));
    }

  });
  std::string linkLines = linkStream.str();
  out_stream << linkLines;
}

void Graph::GfaSerializer::DumpPathFlow(const PathNodeIds& path_flow, usize pathNum, std::string& windowId, std::ostream& out_stream) {
  std::vector<std::string> windowSplit = absl::StrSplit(windowId,':');
  std::string chrom = "chr" + windowSplit[0];
  std::vector<std::string> posSplit = absl::StrSplit(windowSplit[1],'-');
  std::string start = posSplit[0];
  std::string end = posSplit[1];

  std::ostringstream walkStream;
  walkStream << ">" << path_flow.cbegin()->srcId;
  std::for_each(path_flow.cbegin(), path_flow.cend(), [&](const EdgeNodeIds& ep) {
    walkStream << absl::StreamFormat(">%d",ep.dstId);
    //out_stream << absl::StreamFormat("%d -> %d [color=\"%f 0.85 0.85\"]", ep.srcId, ep.dstId, hue);
  });
  out_stream << absl::StreamFormat("W\tsampleName??\t%d\t%s\t%s\t%s\t",pathNum, chrom, start, end);
  std::string walk = walkStream.str();
  out_stream << walk << "\n";
}

auto Graph::GfaSerializer::OppositeStrandSequence(absl::string_view seq) -> std::string {
  std::string result(seq.cbegin(), seq.cend());
  std::transform(result.begin(), result.end(), result.begin(), [](const char& b) { return utils::RevComp(b); });
  return result;
}
}  // namespace lancet2
