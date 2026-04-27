# VarBridge: Fast lifting of variants from a source to a target genome using a source-to-target genome alignment

VarBridge lifts variants from a source to a target genome using an alignment of the source to the target genome.

## Building from source

VarBridge can be built from source using a recursive clone and make. VarBridge depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/VarBridge.git`

`cd VarBridge/`

`make all`

## Usage

For a new telomere-to-telomere assembly (asm) and variants of a sample S1 against that assembly (asm.bcf) the lifting to GRCh38 requires an alignment of asm against GRCh38 (asm_to_hg38.bam), e.g.:

`minimap2 -a -x asm5 --cs hg38.fa asm.fa | samtools sort -O bam -o asm_to_hg38.bam`

With such an alignment, VarBridge lifts variants using:

`varbridge lift -g hg38.fa -s S1 -a asm.bcf asm_to_hg38.bam`

You can redirect the output to a VCF file and also output all non-liftable variants with their closest interval in GRCh38.

`varbridge lift -g hg38.fa -b out.bed -o out.vcf -s S1 -a asm.bcf asm_to_hg38.bam`

The output VCF is unsorted and thus, you need to sort it using bcftools for indexing.

`bcftools sort out.vcf > sorted.vcf`

Variants where the ALT allele in the assembly is the REF allele in the target genome (GRCh38) get the flag REF_ALT_SWAP. Please note that homozygous alternative variants on the assembly with a 1/1 genotype are then converted to 0/0 for the target genome if there is a REF_ALT_SWAP. You can use bcftools to filter such non-variant sites.

`bcftools view --min-ac 1 sorted.vcf` 


## License

VarBridge is free and open source (BSD). Consult the accompanying [LICENSE](https://github.com/tobiasrausch/VarBridge/blob/master/LICENSE) file for more details.


## Credits

[HTSlib](https://github.com/samtools/htslib) is heavily used for alignment processing. [Boost](https://www.boost.org/) for various data structures and algorithms and [Edlib](https://github.com/Martinsos/edlib) for alignments of the REF and ALT allele to the target genome.

