#ifndef LIFT_H
#define LIFT_H

#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h>
#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include "edlib.h"
#include "variants.h"
#include "util.h"

namespace varbridge {


  struct LiftConfig {
    typedef std::map<std::string, uint32_t> TChromMap;
    
    bool multiLift;
    uint16_t minMapQual;
    int32_t win;
    int32_t extra;
    std::string sample;
    boost::filesystem::path vcffile;
    boost::filesystem::path outfile;
    boost::filesystem::path bedfile;
    boost::filesystem::path bamfile;
    boost::filesystem::path genome;
    TChromMap vcfMap;
  };

  struct AlignSegment {
    int32_t asm_start;
    int32_t asm_end;
    int32_t hg38_chr;
    int32_t hg38_min;
    int32_t hg38_max;
  };


  // Map an assembly position to a target genome position.
  inline int32_t
  queryToRefPos(uint32_t const* cigar, int32_t n_cigar, int32_t refStart, int32_t queryPos) {
    int32_t rp = refStart;
    int32_t qp = 0;
    for (int32_t i = 0; i < n_cigar; ++i) {
      int op = bam_cigar_op(cigar[i]);
      int32_t len = bam_cigar_oplen(cigar[i]);
      if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
	if (queryPos < qp + len) return rp + (queryPos - qp);
	qp += len; rp += len;
      } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
	if (queryPos < qp + len) return rp;
	qp += len;
      } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
	rp += len;
      }
    }
    return rp;
  }


  // Left-normalize a variant allele pair against the reference sequence.
  inline void
  leftNormalize(char const* refseq, int32_t /*chromLen*/, int32_t& pos, std::string& ref, std::string& alt) {
    // Trim common tail
    while ((ref.size() > 1) && (alt.size() > 1) && (ref.back() == alt.back())) {
      ref.pop_back();
      alt.pop_back();
    }
    // Prepend if empty
    while ( ( (ref.empty()) || (alt.empty()) ) && (pos > 0) ) {
      --pos;
      char c = (char)std::toupper((unsigned char)refseq[pos]);
      ref.insert(ref.begin(), c);
      alt.insert(alt.begin(), c);
    }
    // Trim common prefix
    while ((ref.size() > 1) && (alt.size() > 1) && (ref.front() == alt.front())) {
      ref.erase(ref.begin());
      alt.erase(alt.begin());
      ++pos;
    }
    // Left-shift: while the rightmost base of both alleles equals the base just left of pos
    while ((pos > 0) && (!ref.empty()) && (!alt.empty()) && (ref.back() == alt.back())) {
      --pos;
      char c = (char)std::toupper((unsigned char)refseq[pos]);
      ref.pop_back(); alt.pop_back();
      ref.insert(ref.begin(), c);
      alt.insert(alt.begin(), c);
    }
  }

  template<typename TConfig, typename TVariant>
  inline int32_t
  liftVariants(TConfig& c, std::vector<TVariant>& variants) {
    typedef std::vector<TVariant> TVariantVec __attribute__((unused));

    // Open output stream
    std::ofstream outf;
    std::ostream* out;
    if (c.outfile.string() == "-") out = &std::cout;
    else {
      outf.open(c.outfile.string());
      if (!outf.is_open()) {
	std::cerr << "Cannot open output file: " << c.outfile.string() << std::endl;
	return 1;
      }
      out = &outf;
    }

    // Open BAM
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // VCF header
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::gregorian::date today = now.date();
    *out << "##fileformat=VCFv4.2\n";
    *out << "##fileDate=" << boost::gregorian::to_iso_string(today) << "\n";
    *out << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    *out << "##INFO=<ID=LIFT_SRC,Number=1,Type=String,Description=\"Source variant in assembly: contig:pos:ref:alt\">\n";
    *out << "##INFO=<ID=EDLIB_EDIST,Number=1,Type=Integer,Description=\"Edit distance of liftover alignment window\">\n";
    *out << "##INFO=<ID=REF_ALT_SWAP,Number=0,Type=Flag,Description=\"ALT assembly allele is REF allele in target genome\">\n";
    *out << "##INFO=<ID=REVERSE,Number=0,Type=Flag,Description=\"Assembly contig of this variant aligns in reverse to target genome\">\n";
    *out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    for (int32_t i = 0; i < hdr->n_targets; ++i) *out << "##contig=<ID=" << hdr->target_name[i] << ",length=" << hdr->target_len[i] << ">\n";
    *out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << c.sample << "\n";

    // Track lift statistics
    std::vector<bool> ever_lifted(variants.size(), false);
    std::map<int32_t, std::vector<AlignSegment>> alignSegments; // Alignment segments per contig

    // Parse BAM alignments
    int32_t refIndex = -1;
    char* seq = NULL;
    int32_t chromLen = 0;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
      if ((rec->core.qual < c.minMapQual) || (rec->core.tid < 0)) continue;

      // Any variants on this contig?
      std::string ctgname = bam_get_qname(rec);
      auto vcfIt = c.vcfMap.find(ctgname);
      if (vcfIt == c.vcfMap.end()) continue;
      int32_t rid = (int32_t)vcfIt->second;

      // Leading hard-clip
      uint32_t* cigar = bam_get_cigar(rec);
      int32_t hc_start = 0;
      if ((rec->core.n_cigar > 0) && (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP)) hc_start = (int32_t)bam_cigar_oplen(cigar[0]);

      // Cigar parsing
      int32_t gp = rec->core.pos;
      int32_t sp = 0;
      int32_t seqStart = -1;
      int32_t seqEnd = -1;
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	int op = bam_cigar_op(cigar[i]);
	int32_t len = bam_cigar_oplen(cigar[i]);
	if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
	  if (seqStart == -1) seqStart = sp;
	  gp += len; sp += len; seqEnd = sp;
	} else if (op == BAM_CINS) {
	  if (seqStart == -1) seqStart = sp;
	  sp += len; seqEnd = sp;
	} else if (op == BAM_CDEL) {
	  if (seqStart == -1) seqStart = sp;
	  gp += len;
	} else if (op == BAM_CSOFT_CLIP) {
	  sp += len;
	} else if (op == BAM_CREF_SKIP) {
	  gp += len;
	} else if (op == BAM_CHARD_CLIP) {
	  sp += len;
	} else {
	  std::cerr << "Unknown CIGAR op " << op << std::endl;
	  return 1;
	}
      }
      if (seqStart == -1) continue;

      // Reverse alignment?
      bool fwd = !(rec->core.flag & BAM_FREVERSE);
      if (!fwd) {
	int32_t tmp = seqStart;
	seqStart = sp - seqEnd;
	seqEnd = sp - tmp;
      }

      // Store alignment segments for non-liftable variants.
      {
	int32_t sp_lo = fwd ? seqStart : (sp - seqEnd);
	int32_t sp_hi = fwd ? (seqEnd - 1) : (sp - seqStart - 1);
	int32_t hg38_lo = queryToRefPos(cigar, (int32_t)rec->core.n_cigar, rec->core.pos, sp_lo);
	int32_t hg38_hi = queryToRefPos(cigar, (int32_t)rec->core.n_cigar, rec->core.pos, sp_hi);
	AlignSegment seg;
	seg.asm_start = seqStart;
	seg.asm_end   = seqEnd;
	seg.hg38_chr  = rec->core.tid;
	seg.hg38_min  = std::min(hg38_lo, hg38_hi);
	seg.hg38_max  = std::max(hg38_lo, hg38_hi) + 1;
	alignSegments[rid].push_back(seg);
      }

      if (rec->core.l_qseq <= 0) continue; // No sequence, skip

      // Lazy loading of target genome
      if (rec->core.tid != refIndex) {
	if ((refIndex != -1) && (seq != NULL)) { free(seq); seq = NULL; }
	refIndex = rec->core.tid;
	int32_t seqlen = -1;
	std::string tname(hdr->target_name[refIndex]);
	chromLen = (int32_t)hdr->target_len[refIndex];
	seq = faidx_fetch_seq(fai, tname.c_str(), 0, chromLen, &seqlen);
	if (seq == NULL) {
	  std::cerr << "Failed to fetch sequence for " << tname << std::endl;
	  refIndex = -1; continue;
	}
      }

      // Get sequence
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

      // Iterate over variants within the aligned span of this contig
      auto itVar = std::lower_bound(variants.begin(), variants.end(), TVariant(rid, seqStart));
      while ((itVar != variants.end()) && (itVar->chr == rid) && (itVar->pos < seqEnd)) {
	std::size_t varIdx = (std::size_t)(itVar - variants.begin());
	int32_t varPos = itVar->pos;
	int32_t ref_len = (int32_t) itVar->ref.size();

	// Lift variants once or for all alignments?
	if (( (ever_lifted[varIdx]) && (!c.multiLift) ) || (varPos < seqStart) || (varPos + ref_len > seqEnd)) {
	  ++itVar; continue;
	}

	// Compute where the REF allele starts in the target genome
	int32_t var_start = fwd ? (varPos - hc_start) : (sp - varPos - ref_len - hc_start);
	if ( (var_start < 0) || (var_start + ref_len > (int32_t) sequence.size() ) ) {
	  ++itVar; continue;
	}

	// Assembly window centred on the REF allele (clamped to stored sequence)
	int32_t asm_s = std::max(0, var_start - c.win);
	int32_t asm_e = std::min((int32_t) sequence.size(), var_start + ref_len + c.win);
	if (asm_s >= asm_e) { ++itVar; continue; }
	std::string asm_query = sequence.substr(asm_s, asm_e - asm_s);
	std::transform(asm_query.begin(), asm_query.end(), asm_query.begin(), ::toupper);
	int32_t var_in_query = var_start - asm_s;  // Offset of REF start in query

	// Approximate GRCh38 position: walk CIGAR in sp-space.
	int32_t queryPos = fwd ? varPos : (sp - varPos - ref_len);
	int32_t approxGP = queryToRefPos(cigar, (int32_t)rec->core.n_cigar, rec->core.pos, queryPos);

	// Target genome window (with extra flank)
	int32_t hg38_s = std::max(0, approxGP - c.win - c.extra);
	int32_t hg38_e = std::min(chromLen, approxGP + ref_len + c.win + c.extra);
	if (hg38_s >= hg38_e) { ++itVar; continue; }
	std::string hg38_target(seq + hg38_s, seq + hg38_e);
	std::transform(hg38_target.begin(), hg38_target.end(), hg38_target.begin(), ::toupper);

	// Align alleles
	EdlibAlignResult res = edlibAlign(asm_query.c_str(), (int)asm_query.size(), hg38_target.c_str(), (int)hg38_target.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
	if ((res.status != EDLIB_STATUS_OK) || (res.alignment == NULL)) {
	  edlibFreeAlignResult(res); ++itVar; continue;
	}

	// Walk the alignment
	int32_t tIdx = (int32_t)infixStart(res);
	int32_t qIdx = 0;
	int32_t t_var = -1;
	for (int32_t ai = 0; ai < res.alignmentLength; ++ai) {
	  if (qIdx == var_in_query) { t_var = tIdx; break; }
	  unsigned char op = res.alignment[ai];
	  if ((op == EDLIB_EDOP_MATCH) || (op == EDLIB_EDOP_MISMATCH)) { ++qIdx; ++tIdx; }
	  else if (op == EDLIB_EDOP_INSERT)  { ++qIdx; }   // base in query only
	  else if (op == EDLIB_EDOP_DELETE)  { ++tIdx; }   // base in target only
	}
	if ((t_var < 0) && (qIdx == var_in_query)) t_var = tIdx; 
	int32_t edist = res.editDistance;
	edlibFreeAlignResult(res);
	if (t_var < 0) { ++itVar; continue; }
	int32_t orig_lifted_gp = hg38_s + t_var;  // 0-based position in GRCh38 (before normalization)
	int32_t lifted_gp = orig_lifted_gp;
	if ((lifted_gp < 0) || (lifted_gp + ref_len > chromLen)) { ++itVar; continue; }

	// Build GRCh38 REF and ALT alleles.
	std::string hg38_ref(seq + lifted_gp, seq + lifted_gp + ref_len);
	std::transform(hg38_ref.begin(), hg38_ref.end(), hg38_ref.begin(), ::toupper);
	std::string hg38_alt = itVar->alt;
	std::transform(hg38_alt.begin(), hg38_alt.end(), hg38_alt.begin(), ::toupper);
	if (!fwd) reverseComplement(hg38_alt);
	// Left-normalize 
	if (hg38_ref != hg38_alt) leftNormalize(seq, chromLen, lifted_gp, hg38_ref, hg38_alt);
	if (hg38_ref == hg38_alt) {
	  // Assembly ALT maps to hg38 REF
	  std::string hg38_ref2(seq + orig_lifted_gp, seq + orig_lifted_gp + ref_len);
	  std::transform(hg38_ref2.begin(), hg38_ref2.end(), hg38_ref2.begin(), ::toupper);
	  std::string hg38_alt2 = itVar->ref;
	  std::transform(hg38_alt2.begin(), hg38_alt2.end(), hg38_alt2.begin(), ::toupper);
	  if (!fwd) reverseComplement(hg38_alt2);
	  int32_t lifted_gp2 = orig_lifted_gp;
	  if (hg38_ref2 != hg38_alt2) leftNormalize(seq, chromLen, lifted_gp2, hg38_ref2, hg38_alt2);
	  // Swap genotypes
	  int32_t new_gtA1 = (itVar->gtA1 == 0) ? 1 : 0;
	  int32_t new_gtA2 = (itVar->gtA2 == 0) ? 1 : 0;
	  *out << hdr->target_name[refIndex] << '\t' << (lifted_gp2 + 1) << '\t' << '.' << '\t' << hg38_ref2 << '\t' << hg38_alt2 << '\t';
	  if (bcf_float_is_missing(itVar->qual)) *out << '.';
	  else *out << itVar->qual;
	  *out << '\t' << itVar->filter << '\t' << "LIFT_SRC=" << ctgname << ':' << (varPos + 1) << ':' << itVar->ref << ':' << itVar->alt << ";EDLIB_EDIST=" << edist << ";REF_ALT_SWAP";
	  if (rec->core.flag & BAM_FREVERSE) *out << ";REVERSE";
	  *out << '\t' << "GT" << '\t' << new_gtA1 << '/' << new_gtA2 << std::endl;
	  ever_lifted[varIdx] = true;
	  ++itVar; continue;
	}

	// VCF record
	*out << hdr->target_name[refIndex] << '\t' << (lifted_gp + 1) << '\t' << '.' << '\t' << hg38_ref << '\t' << hg38_alt << '\t';
	if (bcf_float_is_missing(itVar->qual)) *out << '.';
	else *out << itVar->qual;
	*out << '\t' << itVar->filter << '\t' << "LIFT_SRC=" << ctgname << ':' << (varPos + 1) << ':' << itVar->ref << ':' << itVar->alt << ";EDLIB_EDIST=" << edist;
	if (rec->core.flag & BAM_FREVERSE) *out << ";REVERSE";
	*out << '\t' << "GT" << '\t' << itVar->gtA1 << '/' << itVar->gtA2 << std::endl;
	ever_lifted[varIdx] = true;
	++itVar;
      }
    }

    // Clean up
    bam_destroy1(rec);
    if (seq != NULL) free(seq);
    fai_destroy(fai);

    // BED output for non-liftable variants
    if (!c.bedfile.empty()) {
      std::ofstream bedout(c.bedfile.string());
      if (!bedout.is_open()) {
	std::cerr << "Cannot open BED output file: " << c.bedfile.string() << std::endl;
	bam_hdr_destroy(hdr);
	hts_idx_destroy(idx);
	sam_close(samfile);
	return 1;
      }

      // Build reverse map
      std::map<int32_t, std::string> ridToCtg;
      for (auto const& kv : c.vcfMap) ridToCtg[(int32_t)kv.second] = kv.first;

      // Sort each contig's segments by asm_start for binary search
      for (auto& kv : alignSegments) std::sort(kv.second.begin(), kv.second.end(), [](AlignSegment const& a, AlignSegment const& b) { return a.asm_start < b.asm_start; });

      // Iterate all variants
      for (size_t i = 0; i < variants.size(); ++i) {
	if (ever_lifted[i]) continue;
	auto const& v = variants[i];
	int32_t ref_len = (int32_t)v.ref.size();

	auto segIt = alignSegments.find(v.chr);
	if (segIt == alignSegments.end()) continue;
	auto const& segs = segIt->second;

	// Find right segment
	auto right_it = std::lower_bound(segs.begin(), segs.end(), v.pos + ref_len, [](AlignSegment const& s, int32_t val) { return s.asm_start < val; });
	// Find left segment
	auto left_bound = std::lower_bound(segs.begin(), segs.end(), v.pos, [](AlignSegment const& s, int32_t val) { return s.asm_start < val; });
	// Search backward for a segment whose asm_end <= v.pos
	auto left_it = segs.end();
	for (auto it = left_bound; it != segs.begin(); ) {
	  --it;
	  if (it->asm_end <= v.pos) { left_it = it; break; }
	}
	if ((left_it == segs.end()) || (right_it == segs.end())) continue;
	if (left_it->hg38_chr != right_it->hg38_chr) continue;

	// Output flanking positions
	int32_t bed_s = std::min(left_it->hg38_max, right_it->hg38_max);
	int32_t bed_e = std::max(left_it->hg38_min, right_it->hg38_min);
	if (bed_s >= bed_e) continue; 
	std::string ctgname = ridToCtg.count(v.chr) ? ridToCtg[v.chr] : std::to_string(v.chr);
	bedout << hdr->target_name[left_it->hg38_chr] << '\t' << bed_s << '\t' << bed_e << '\t' << ctgname << ':' << (v.pos + 1) << ':' << v.ref << ':' << v.alt << std::endl;
      }
    }

    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return 0;
  }


  template<typename TConfig>
  inline int32_t
  runLift(TConfig& c) {

#ifdef PROFILE
    ProfilerStart("lift.prof");
#endif

    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Load variants" << std::endl;
    std::vector<BiallelicVariant> variants;
    if (!_loadVariants(c, variants)) {
      std::cerr << "Loading variants failed!" << std::endl;
      return 1;
    } else {
      if (variants.empty()) {
	std::cerr << "No variants found for " << c.sample << ". Please check the VCF sample name!" << std::endl;
	return 1;
      }
    }
    std::sort(variants.begin(), variants.end());

    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Lift variants" << std::endl;;   
    if (liftVariants(c, variants) != 0) {
      std::cerr << "Lifting variants failed!" << std::endl;
      return 1;
    }

#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cerr << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;;
    
    return 0;
  }

  int lift(int argc, char **argv) {
    LiftConfig c;
   
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "BCF sample name")
      ("variants,a", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "BCF input file (variants on assembly)")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "BCF output file (variants on target genome)")
      ;
    
    boost::program_options::options_description disc("Lift options");
    disc.add_options()
      ("map-qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. contig mapping quality")
      ("window,n", boost::program_options::value<int32_t>(&c.win)->default_value(50), "realignment window")
      ("extra,e", boost::program_options::value<int32_t>(&c.extra)->default_value(50), "extra target genome flank")
      ("multi-lift,m", boost::program_options::bool_switch(&c.multiLift)->default_value(false), "lift variant for all covering alignments")
      ("bed-file,b", boost::program_options::value<boost::filesystem::path>(&c.bedfile), "BED output file for non-liftable variants")
      ;

    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(disc).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(disc);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
      std::cerr << std::endl;
      std::cerr << "Usage: varbridge " << argv[0] << " [OPTIONS] -g <hg38.fa> -a <assembly.variants.bcf> <assembly_to_hg38.bam> " << std::endl;
      std::cerr << visible_options << "\n";
      return 0;
    }

    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.genome.string().c_str()) == -1) {
	  std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.genome.string().c_str());
      }
      fai_destroy(fai);
    }
   
    // Check input files
    {
      if (!(boost::filesystem::exists(c.bamfile) && boost::filesystem::is_regular_file(c.bamfile) && boost::filesystem::file_size(c.bamfile))) {
	std::cerr << "Alignment file is missing: " << c.bamfile.string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bamfile.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.bamfile.string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamfile.string() << std::endl;
	return 1;
      }
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (!faidx_has_seq(fai, tname.c_str())) {
	  std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	  return 1;
	}
      }
      fai_destroy(fai);
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
    
    // Check outfile
    if (!vm.count("outfile")) c.outfile = "-";
    else {
      if (c.outfile.string() != "-") {
	if (!_outfileValid(c.outfile)) return 1;
      }
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "varbridge ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    return runLift(c);
  }



}

#endif
