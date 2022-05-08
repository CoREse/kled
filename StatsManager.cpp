/*
 * StatOutputer.cpp
 *
 *  Created on: 2020年5月1日
 *      Author: fenghe
 *	Edited by: CRE 20210705
 */


//
// Manta - Structural Variant and Indel Caller
// Copyright (c) 2013-2019 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

#include "StatsManager.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <unistd.h>

typedef struct R_region{//break point candidate
uint16_t chr_ID;
int  st_pos;
int  ed_pos;
}R_region;

typedef struct Bam_file{
bool _is_record_set;//set to 1 when _brec has data
htsFile* _hfp;//the file of BAM/CRAM
bam_hdr_t* _hdr;//header for BAM/CRAM file
hts_idx_t* _hidx;//index for bam/cram file
hts_itr_t* _hitr;//Iterator for bam/cram file
bam1_t _brec;//current BAM record

// track for debug only:
unsigned _record_no;
char * _stream_name;
bool _is_region;
char* _region;
}Bam_file;

//to = pre + suf; return to
char* strcmb(char* to, char *pre, char * suf)
{
	strcpy(to, pre);
	strcat(to, suf);
	return to;
}

uint8_t bam_map_qual(bam1_t* _bp) 	{	return _bp->core.qual;}

bool fexist(const char *filename)
{
	if(access(filename, 0) == 0)
		return true;
	else
		return false;
}

const uint32_t* bam_raw_cigar(bam1_t* _bp){
	return bam_get_cigar(_bp);
}

unsigned bam_n_cigar(bam1_t* _bp){
	return _bp->core.n_cigar;
}

/// \brief Test if this read contains an 'SA' tag, used to annotate split read alignments
///
/// \return True if the 'SA' tag is found
bool bam_isSASplit(bam1_t* _bp)
{
	static char satag[] = { 'S', 'A' };
	return (NULL != bam_get_string_tag(_bp, satag));
}

/// \brief Test if the read is supplemental, using a more liberal community criteria to define
/// 'supplemental' compared to that from the BAM spec.
///
/// Reads are considered supplemental if either:
/// 1. The 'supplemental' bit is set in the bam record.
/// 2. The 'secondary' bit is set in the bam record and the record contains an 'SA' tag.
///
/// The second condition supports the common workaround typified by bwamem's '-M' option,
/// which allows split reads to be added to the alignment without creating BAM's which could break
/// on older tools.
///
/// \return True if this read is treated as supplemental
bool bam_isNonStrictSupplement(bam1_t* _bp)
{
	if (bam_is_supplementary(_bp))
		return true;
	if (!bam_is_secondary(_bp))
		return false;
	return bam_isSASplit(_bp);
}

/// \brief Test if the read is supplemental, using a more liberal community criteria to define
/// 'supplemental' compared to that from the BAM spec.
///
/// Reads are considered supplemental if either:
/// 1. The 'supplemental' bit is set in the bam record.
/// 2. The 'secondary' bit is set in the bam record and the record contains an 'SA' tag.
///
/// The second condition supports the common workaround typified by bwamem's '-M' option,
/// which allows split reads to be added to the alignment without creating BAM's which could break
/// on older tools.
///
/// \return True if this read is treated as supplemental
bool isNonStrictSupplement(bam1_t* b)
{
  if (bam_is_supplementary(b)) return true;
  if (!bam_is_secondary(b)) return false;
  return bam_isSASplit(b);
}

void bam_load_index(Bam_file * bf)
{
  if (NULL != bf->_hidx)//return directly when header file exist
	  return;

  char *bam_cram_file_name = bf->_hfp->fn;

  // check whether index files exist
  char index_name[1024];
  if(
		  (!fexist( strcmb(index_name, bam_cram_file_name, ".bai"))) &&
		  (!fexist( strcmb(index_name, bam_cram_file_name, ".csi"))) &&
		  (!fexist( strcmb(index_name, bam_cram_file_name, ".crai"))))
  {
	  //build index for cram/bam files
	  fprintf(stderr, "BAM/CRAM index is not available for file %s, now building new index file for it.\n", bam_cram_file_name);
	  sam_index_build3(bam_cram_file_name, NULL, 0, 4);
	  //bam_build_index(bam_cram_file_name);
	  fprintf(stderr, "End for building BAM/CRAM index for file %s\n", bam_cram_file_name);
  }

  bf->_hidx = sam_index_load(bf->_hfp, bam_cram_file_name);
  assert(NULL != bf->_hidx);//, "BAM/CRAM index is not available for file");
}

void resetRegion_ID(Bam_file * bf, R_region *region)
{
  if (NULL != bf->_hitr)
	  hts_itr_destroy(bf->_hitr);
  bam_load_index(bf);

  assert(region->chr_ID >= 0);//, "Invalid region specified for BAM/CRAM file");

  bf->_hitr = sam_itr_queryi(bf->_hidx, region->chr_ID, region->st_pos, region->ed_pos);
  assert(bf->_hitr != NULL);//, "Failed to fetch region specified for BAM/CRAM file");
  bf->_is_region = true;
  //bf->_region.clear();//todo::

  bf->_is_record_set = false;
  bf->_record_no     = 0;
}

void bam_file_close(Bam_file * bf){
	hts_close(bf->_hfp);
	bam_hdr_destroy(bf->_hdr);
	hts_idx_destroy(bf->_hidx);
	hts_itr_destroy(bf->_hitr);
	if(bf->_stream_name) free(bf->_stream_name);
	if(bf->_region) free(bf->_region);
	memset(bf, 0, sizeof(Bam_file));
}

//read a record and stored in "_brec"
bool bam_next(Bam_file * bf)
{
	if (NULL == bf->_hfp)
		return false;
	int ret;
	if (NULL == bf->_hitr)
	{
		ret = sam_read1(bf->_hfp, bf->_hdr, &(bf->_brec));
		// Semi-documented sam_read1 API: -1 is expected read failure at end of stream, any other negative value
		assert(ret >= -1);//, "Unexpected return value from htslib sam_read1 function while attempting to read BAM/CRAM file:\n");
	}
	else
	{
		ret = sam_itr_next(bf->_hfp, bf->_hitr, &(bf->_brec));
		assert(ret >= -1);//, "Unexpected return value from htslib sam_read1 function while attempting to read BAM/CRAM file:\n");
	}
	bf->_is_record_set = (ret >= 0);
	if (bf->_is_record_set) bf->_record_no++;

	return bf->_is_record_set;
}

void parse_bam_region(char* region, char* chrom, int32_t* begin_pos, int32_t* end_pos)
{
	char* chrom_;
	int32_t begin_pos_ = 0;
	int32_t end_pos_ = 0;

	// make first split:
	char region_cp[512];
	strcpy(region_cp, region);
	char* afterChrom = strchr(region, ':');
    if (NULL != afterChrom)
    {
    	chrom_ = region;
    	afterChrom[0] = '\0';
    	assert(afterChrom[1] != '\0');//, "");
    	afterChrom++;
    }

    bool isWholeChrom = (NULL == afterChrom);
    if (!isWholeChrom)
    {
    	// make second split
    	char *tokens = NULL;
    	tokens = strtok(afterChrom, "-") - 1;
    	begin_pos_ = strtoul(tokens,NULL, 10);
    	tokens = strtok(NULL, "\0");
    	end_pos_ = strtoul(tokens,NULL, 10);
    	//"Can't parse begin and end positions from bam_region"
    }
    if(begin_pos_ < 0 || begin_pos_ > end_pos_)
    	isWholeChrom = true;
    assert(chrom_ != NULL);//, "Can't parse contig name from bam_region ");
	chrom     = chrom_;
    if (isWholeChrom)
    {
    	*begin_pos = 0;
    	*end_pos   = 0xffffffff;
    }
    else
    {
    	*begin_pos = begin_pos_;
    	*end_pos   = end_pos_;
    }
}

#define MIN(a,b) (((a) < (b))?(a):(b))

void region_string2R_region(
    bam_hdr_t* header,
	char* region_str,
	R_region *region)
{
  assert(NULL != header);//, "");
  assert(NULL != region_str);//, "");
  char * chrom = NULL;
  parse_bam_region(region_str, chrom, &(region->st_pos), &(region->ed_pos));

  region->chr_ID = bam_name2id(header, chrom);
  assert(region->chr_ID >= 0);//, "Contig [%s] from bam_region [%s] not found in BAM/CRAM header",chrom, region_str, NULL);
  region->ed_pos = MIN(region->ed_pos, header->target_len[region->chr_ID]);
}

void resetRegion_char(Bam_file * bf, char* region_str)
{
	R_region region;
	region_string2R_region(bf->_hdr, region_str, &region);
    resetRegion_ID(bf,  &region);
    bf->_region = region_str;
}

void bam_file_open(const char* filename, const char* referenceFilename, char* region, Bam_file * bf)
{
	assert(filename != NULL);//, "Can't initialize bam_streamer with empty filename\n");
	assert(*filename != '\0');//, "Can't initialize bam_streamer with empty filename\n");
	//try open:
	FILE* try_open = fopen(filename, "r");
	fclose(try_open);
	memset(bf, 0, sizeof(Bam_file));
	bf->_hfp = hts_open(filename, "rb");
	assert(bf->_hfp != NULL);//, "Failed to open SAM/BAM/CRAM file for reading\n");

	//set reference file
	if (NULL != referenceFilename)
	{
		char referenceFilenameIndex[128];
		strcpy(referenceFilenameIndex, referenceFilename);
		strcat(referenceFilenameIndex, ".fai");
		int ret = hts_set_fai_filename(bf->_hfp, referenceFilenameIndex);
		assert(ret == 0);//, "Failed to use reference for BAM/CRAM file");
	}
	bf->_hdr = sam_hdr_read(bf->_hfp);
	assert(bf->_hdr != NULL);//, "Failed to parse header from SAM/BAM/CRAM file");

	//set region
	if (NULL == region)
	{
		// setup to read the whole BAM file by default if resetRegion() is not called:
		if (bf->_hdr->n_targets)
			// parse any contig name so that header->hash is created; ignore returned tid value, so doesn't matter if fake name exists
		bam_name2id(bf->_hdr, "fake_name");
	}
	else
		// read a specific region of the bam file:
		resetRegion_char(bf, region);
}

inline bool is_segment_align_match(const int id)
{
  switch (id) {
  case CIGAR_MATCH:
  case CIGAR_SEQ_MATCH:
  case CIGAR_SEQ_MISMATCH:
    return true;
  default:
    return false;
  }
}

bool readFilteredAlignment(bam1_t* b)
{
	bool isMatched = false;
	bool isSkip = false;
	bool isClipped = false;
	bool reverse = !bam_is_fwd_strand(b);
	const uint32_t* bam_cigar = bam_raw_cigar(b); unsigned n_cigar = bam_n_cigar(b);
	for (int i = 0; i < n_cigar; ++i)
	{
		uint32_t cigar = reverse?bam_cigar[n_cigar - i - 1]:bam_cigar[i];
		int len = (cigar >> BAM_CIGAR_SHIFT);
		int type = (int)(1 + (cigar & BAM_CIGAR_MASK));

		if (is_segment_align_match(type)) {
			if (isClipped) return true;
			isMatched = true;
		}
		else if (type == CIGAR_SKIP) {
			if (isSkip) return true;
			isSkip = true;
		}
		else if (type == CIGAR_SOFT_CLIP)
			isClipped = true;
		else
			return true;
	}
	return (!isMatched);
}

bool bam_is_mapped_pair(bam1_t* b)
{
  if (!bam_is_paired(b)) return false;
  if (bam_is_unmapped(b) || bam_is_mate_unmapped(b)) return false;
  return true;
}

bool hasRefSkip(bam1_t* b)
{
	const uint32_t* bam_cigar = bam_raw_cigar(b); unsigned n_cigar = bam_n_cigar(b);
	for (int i = 0; i < n_cigar; ++i)
		 if ((int)(1 + (bam_cigar[i] & BAM_CIGAR_MASK)) == CIGAR_SKIP)
			 return true;
	return false;
}

bool isReadFilteredCore(bam1_t* b)
{
  if (bam_is_filter(b))
    return true;
  else if (bam_is_dup(b))
    return true;
  // supplementary reads without SA tag
  else if (bam_is_supplementary(b) && (!bam_isSASplit(b)))
    return true;
  else
  {
    // hack to work with bwamem '-M' formatting,
    // keep secondary reads when they contain an SA tag
    if (bam_is_secondary(b))
    {
      if (!bam_isSASplit(b)) return true;
    }
  }
  return false;
}

char* bam_get_string_tag(bam1_t* _bp, char* tag)
{
  // retrieve the BAM tag
  uint8_t* pTag = bam_aux_get(_bp, tag);
  if (!pTag) return NULL;

  // skip tags that are not encoded as a null-terminated string
  if (pTag[0] != 'Z') return NULL;
  ++pTag;

  return (char*)pTag;
}

bool ReadPairDepthFilter::isFilterRead(bam1_t *b) {
		static const unsigned maxPosCount(1);

		if (b->core.tid != _lastTargetId) {
			_goodMates.clear();
			_lastTargetId = b->core.tid;
			_posCount = 0;
			_lastPos = b->core.pos;
		} else if (b->core.pos != _lastPos) {
			_posCount = 0;
			_lastPos = b->core.pos;
		}

		// Assert only two reads per fragment
		const unsigned readNum(bam_is_first(b) ? 1 : 2);
		assert(bam_is_second(b) == (readNum == 2));

		// Filter pairs with templateSize 0 (unknown)
		if (b->core.isize == 0)
			return true;

		// sample each read pair once by sampling stats from
		// downstream read only, or whichever read is encountered
		// second if the read and its mate start at the same position:
		const bool isDownstream(b->core.pos > b->core.mpos);
		const bool isSamePos(b->core.pos == b->core.mpos);

		if (isDownstream || isSamePos) {
			const int mateReadNo(bam_is_first(b) ? 2 : 1); //当read是first的时候，mate的No就是2；当read的No是2的时候，mate的No就是1
			mateMap_t::iterator i(_goodMates.find(getKey((char*) bam_get_qname(b), mateReadNo)));

			if (i == _goodMates.end()) {    //i == the last one: item not found
				if (isDownstream)
					return true;
			} else {
				_goodMates.erase(i);
				return false;
			}
		}

		// to prevent high-depth pileups from overly biasing the
		// read stats, we only take maxPosCount read pairs from each start
		// pos. by not inserting a key in goodMates, we also filter
		// the downstream mate:
		if (_posCount >= maxPosCount)
			return true;
		++_posCount;

		// crude mechanism to manage total set memory
		if (_goodMates.size() > maxMateSetSize)
			_goodMates.clear();

		// Ignore pairs where the upstream mate has a refskip, since we cannot
		// compute the correct insert size later when looking at the downstream mate
		// (Or we would have to save the total refskip length here)
		if (hasRefSkip(b))
			return true;

		_goodMates.insert(getKey(b));
		return true;
	}

bool CoreReadFilter::isFilterRead(bam1_t *b) {
	// filter common categories of undesirable reads:
	if (isReadFilteredCore(b))
		return true;
	if (isNonStrictSupplement(b))
		return true;
	if (!bam_is_mapped_pair(b))
		return true;
	if (bam_map_qual(b) == 0)
		return true;
	if (bam_isSASplit(b))
		return true;    // filter any split reads with an SA tag:
	if (readFilteredAlignment(b))
		return true; // remove alignments other than {X}M({Z}N{X2}M)?({Y}S)? (or reverse for reverse strand)
	if (pairFilter.isFilterRead(b))
		return true;  // filter out upstream reads and high depth regions:
	return false;
}


/// this struct exists for the sole purpose of xml output:
struct ReadGroupStatsExporter {

	std::string bamFile;
	std::string readGroup;
	UniqueStats groupStats;
};

void StatsManager::merge(const StatsManager &rhs) {
	const unsigned numGroups(rhs.size());
	for (unsigned i(0); i < numGroups; ++i) {
		const StatLabel &mkey(rhs.getKey(i));
		if (_group.test_key(mkey)) {
			std::cerr << "Can't merge stats set objects with repeated key: '"
					<< mkey << "'";
			exit(EXIT_FAILURE);
		}

		setStats(mkey, rhs.getStats(i));
	}
}

#ifdef READ_GROUPS
	#define USE_RG true
#else
	#define USE_RG false
#endif

void StatsManager::handleBamCramStats(const char *alignmentFilename)
{
	TrackerManager rgManager(USE_RG, alignmentFilename,	defaultStatsFilename);
	Bam_file bf = {0}; bam_file_open(alignmentFilename, referenceFilename, NULL, &bf);
	const bam_hdr_t &header(*(bf._hdr));

	std::vector<ChromInfo> chromList;
	for (int32_t i(0); i < header.n_targets; ++i)
		chromList.push_back(ChromInfo(i, header.target_len[i]));

	bool isStopEstimation(false);
	bool isActiveChrom(true);

	while (isActiveChrom && (!isStopEstimation)) {
		isActiveChrom = false;

		for (auto &chrom : chromList) {
			if (isStopEstimation)
				break; // keep sampling until either the chromosome has been exhuasted or the current chunk has been sufficiently sampled
			bool isFinishedSlice(false);
			R_region region;
			region.chr_ID = chrom.chr_ID;
			while (!isFinishedSlice) {
				const int32_t startPos(chrom.highestPos + 1);
				if (startPos >= chrom.size)
					break;
				region.st_pos = startPos;
				region.ed_pos = chrom.size;
				resetRegion_ID(&bf, &region);

				while (bam_next(&bf)) {
					bam1_t *b = &(bf._brec);
					if (b->core.pos < startPos)
						continue;
					chrom.highestPos = b->core.pos;
					isActiveChrom = true;

					StatsTracker &rgInfo(rgManager.getTracker(b));
					rgInfo.handleReadRecordBasic(b);
					if (rgManager.isFiltered(b))
						continue;
					RGT_RETURN r = rgInfo.handleReadRecordCheck(b);
					if (r == RGT_CONTINUE)
						continue;
					else if (r == RGT_BREAK) {
						chrom.highestPos += std::max(1, chrom.size / 100);
						break;
					} //else do nothing

					isFinishedSlice = rgManager.isFinishedSlice();
					if (!isFinishedSlice)
						continue;
					isStopEstimation = rgManager.isStopEstimation();
					// break from reading the current chromosome
					break;
				}
				// move to next region if no read falling in the current region
				if (chrom.highestPos <= startPos) {
					chrom.highestPos += std::max(1, chrom.size / 100);
				}
			}
		}
	}

	for (auto &val : rgManager.getMap())
		setStats(val.first, val.second.getStats());
	bam_file_close(&bf);
}
