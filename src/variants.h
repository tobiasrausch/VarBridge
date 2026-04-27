#ifndef VARIANTS_H
#define VARIANTS_H

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>


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

  template<typename TConfig, typename TVariant>
  inline bool
  _loadVariants(TConfig& c, std::vector<TVariant>& pV) {
    // Load BCF file
    htsFile* ifile = bcf_open(c.vcffile.c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.c_str());
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
      
      // Only bi-allelic variants
      if (rec->n_allele == 2) {
	bcf_unpack(rec, BCF_UN_ALL);
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
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
    
    return true;
  }

}

#endif
