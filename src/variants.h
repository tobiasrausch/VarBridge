#ifndef VARIANTS_H
#define VARIANTS_H

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>
#include <htslib/vcf.h>


namespace varbridge
{

  struct BiallelicVariant {
    int32_t chr;
    int32_t pos;
    std::string ref;
    std::string alt;
    float qual;
    std::string filter;
    int32_t gtA1;
    int32_t gtA2;

    BiallelicVariant(int32_t const c, int32_t const p) : chr(c), pos(p), ref(""), alt(""), qual(0), filter("."), gtA1(0), gtA2(0) {}
    BiallelicVariant(int32_t const c, int32_t const p, std::string const& r, std::string const& a) : chr(c), pos(p), ref(r), alt(a), qual(0), filter("."), gtA1(0), gtA2(0) {}
    BiallelicVariant(int32_t const c, int32_t const p, std::string const& r, std::string const& a, float const q, std::string const& f) : chr(c), pos(p), ref(r), alt(a), qual(q), filter(f), gtA1(0), gtA2(0) {}
    BiallelicVariant(int32_t const c, int32_t const p, std::string const& r, std::string const& a, float const q, std::string const& f, int32_t const a1, int32_t const a2) : chr(c), pos(p), ref(r), alt(a), qual(q), filter(f), gtA1(a1), gtA2(a2) {}

    bool operator<(BiallelicVariant const& s2) const {
      return ((chr < s2.chr) || ((chr == s2.chr) && (pos < s2.pos)));
    }
  };

  // Structural variant with symbolic ALTs
  struct SvVariant {
    int32_t chr;
    int32_t pos;
    int32_t chr2;
    int32_t svEnd;
    int32_t svLen;
    std::string svType;
    std::string ref;
    std::string alt;
    float qual;
    std::string filter;
    int32_t gtA1;
    int32_t gtA2;
    std::string id;
    std::string chr2Name;
    std::string infoStr;

    SvVariant(int32_t c, int32_t p) : chr(c), pos(p), chr2(c), svEnd(p), svLen(0),
      svType(""), ref(""), alt(""), qual(0), filter("."), gtA1(0), gtA2(0), id(".") {}

    bool operator<(SvVariant const& s2) const {
      return (chr < s2.chr) || ((chr == s2.chr) && (pos < s2.pos));
    }
  };

  template<typename TConfig, typename TVariant>
  inline bool
  _loadVariants(TConfig& c, std::vector<TVariant>& pV) {
    // Load BCF file
    htsFile* ifile = bcf_open(c.vcffile.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    // Find sample
    int32_t sampleIndex = -1;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
      if (hdr->samples[i] == c.sample) sampleIndex = i;
    if (sampleIndex < 0) return false;
        
    // Genotypes
    int ngt = 0;
    int32_t* gt = NULL;

    // Collect bi-allelic variants for all chromosomes
    bcf1_t* rec = bcf_init1();
    int32_t lastId = -1;
    int32_t lastpos = -1;
    while(bcf_read(ifile, hdr, rec) == 0) {
      // Update VCF chromosome map
      if (rec->rid != lastId) {
	lastId = rec->rid;
	std::string chrName = bcf_hdr_id2name(hdr, rec->rid);
	c.vcfMap[chrName] = rec->rid;
      }
      
      // Only bi-allelic variants with exact (non-symbolic) ALT alleles
      if (rec->n_allele == 2) {
	bcf_unpack(rec, BCF_UN_ALL);
	std::string altAllele(rec->d.allele[1]);
	if ((altAllele.front() == '<' && altAllele.back() == '>') || (altAllele.find('[') != std::string::npos) || (altAllele.find(']') != std::string::npos)) continue;
	bcf_get_genotypes(hdr, rec, &gt, &ngt);
	if ((bcf_gt_allele(gt[sampleIndex*2]) != -1) && (bcf_gt_allele(gt[sampleIndex*2 + 1]) != -1) && (!bcf_gt_is_missing(gt[sampleIndex*2])) && (!bcf_gt_is_missing(gt[sampleIndex*2 + 1]))) {
	  int gt_type = bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2 + 1]);
	  if (gt_type > 0) {
	    if (rec->pos != lastpos) {
	      std::string ref = std::string(rec->d.allele[0]);
	      std::string alt = std::string(rec->d.allele[1]);
	      std::string filterStr = ".";
	      if (rec->d.n_flt > 0) {
		filterStr.clear();
		for (int fi = 0; fi < rec->d.n_flt; fi++) {
		  if (fi > 0) filterStr += ";";
		  filterStr += hdr->id[BCF_DT_ID][rec->d.flt[fi]].key;
		}
	      }
	      pV.push_back(TVariant(rec->rid, rec->pos, ref, alt, rec->qual, filterStr, bcf_gt_allele(gt[sampleIndex*2]), bcf_gt_allele(gt[sampleIndex*2 + 1])));
	      lastpos = rec->pos;
	    }
	  }
	}
      }
    }
    bcf_destroy(rec);
    if (gt != NULL) free(gt);
    
    // Close BCF
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    
    return true;
  }

  template<typename TConfig>
  inline bool
  _loadSvVariants(TConfig& c, std::vector<SvVariant>& svV) {
    htsFile* ifile = bcf_open(c.vcffile.c_str(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(ifile);

    // VCF map
    for (int32_t i = 0; i < hdr->n[BCF_DT_CTG]; ++i)
      c.vcfMap[hdr->id[BCF_DT_CTG][i].key] = (uint32_t)i;

    // Find sample to lift
    int32_t sampleIndex = -1;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
      if (hdr->samples[i] == c.sample) sampleIndex = i;
    if (sampleIndex < 0) { bcf_hdr_destroy(hdr); bcf_close(ifile); return false; }

    int ngt = 0;
    int32_t* gt   = NULL;
    int32_t* ivals = NULL;
    int nivals = 0;
    char* svals = NULL;
    int nsvals = 0;
    kstring_t ks = {0, 0, NULL};

    bcf1_t* rec = bcf_init1();
    while (bcf_read(ifile, hdr, rec) == 0) {
      if (rec->n_allele != 2) continue;
      bcf_unpack(rec, BCF_UN_ALL);

      std::string alt(rec->d.allele[1]);
      bool isSymbolic = ((alt.front() == '<') && (alt.back() == '>'));
      bool isBnd = ((alt.find('[') != std::string::npos) || (alt.find(']') != std::string::npos));
      if (!isSymbolic && !isBnd) continue;

      // Require a non-ref genotype
      bcf_get_genotypes(hdr, rec, &gt, &ngt);
      if ((bcf_gt_allele(gt[sampleIndex*2]) == -1) || (bcf_gt_allele(gt[sampleIndex*2+1]) == -1) || bcf_gt_is_missing(gt[sampleIndex*2]) || bcf_gt_is_missing(gt[sampleIndex*2+1])) continue;
      if (bcf_gt_allele(gt[sampleIndex*2]) + bcf_gt_allele(gt[sampleIndex*2+1]) == 0) continue;

      // SVTYPE (required)
      nsvals = 0;
      if (bcf_get_info_string(hdr, rec, "SVTYPE", &svals, &nsvals) <= 0 || svals == NULL) continue;
      std::string svType(svals);

      // END
      int32_t svEnd = rec->pos;
      nivals = 0;
      if (bcf_get_info_int32(hdr, rec, "END", &ivals, &nivals) > 0 && nivals > 0) svEnd = ivals[0] - 1;

      SvVariant sv(rec->rid, rec->pos);
      sv.svEnd = svEnd;
      sv.svType = svType;
      sv.ref = std::string(rec->d.allele[0]);
      sv.alt = alt;
      sv.qual = rec->qual;
      sv.id = std::string(rec->d.id);
      sv.chr2 = rec->rid;
      sv.gtA1 = bcf_gt_allele(gt[sampleIndex*2]);
      sv.gtA2 = bcf_gt_allele(gt[sampleIndex*2+1]);

      if (rec->d.n_flt > 0) {
	sv.filter.clear();
	for (int fi = 0; fi < rec->d.n_flt; fi++) {
	  if (fi > 0) sv.filter += ";";
	  sv.filter += hdr->id[BCF_DT_ID][rec->d.flt[fi]].key;
	}
      }

      // BND: CHR2 and POS2
      if (svType == "BND") {
	nsvals = 0;
	if (bcf_get_info_string(hdr, rec, "CHR2", &svals, &nsvals) > 0 && svals != NULL) {
	  sv.chr2Name = std::string(svals);
	  auto it = c.vcfMap.find(sv.chr2Name);
	  sv.chr2 = (it != c.vcfMap.end()) ? (int32_t)it->second : -1;
	}
	nivals = 0;
	if (bcf_get_info_int32(hdr, rec, "POS2", &ivals, &nivals) > 0 && nivals > 0) sv.svEnd = ivals[0] - 1;
      }

      // SVLEN
      nivals = 0;
      if (bcf_get_info_int32(hdr, rec, "SVLEN", &ivals, &nivals) > 0 && nivals > 0) sv.svLen = std::abs(ivals[0]);
      else sv.svLen = std::max(0, sv.svEnd - sv.pos);

      // Full original INFO string
      ks.l = 0;
      vcf_format(hdr, rec, &ks);
      {
	const char* p = ks.s;
	int tabs = 0;
	while (*p && tabs < 7) { if (*p == '\t') ++tabs; ++p; }
	const char* q = p;
	while (*q && *q != '\t' && *q != '\n') ++q;
	sv.infoStr = std::string(p, q);
      }
      svV.push_back(sv);
    }

    bcf_destroy(rec);
    if (gt    != NULL) free(gt);
    if (ivals != NULL) free(ivals);
    if (svals != NULL) free(svals);
    if (ks.s  != NULL) free(ks.s);
    bcf_hdr_destroy(hdr);
    bcf_close(ifile);
    return true;
  }

}

#endif
