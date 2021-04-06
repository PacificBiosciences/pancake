/**
 * Copyright (c) 2021, Huxelerate S.r.l.
 * All rights reserved.
 *
 * Provides a set of APIs to interact with the HUGenomic ksw2_score core
 * that performs accelerated sequence alignment on FPGA.
 *
 * The core functionality is based on the ksw2_extd method from the ksw2 library: 
 * https://github.com/lh3/ksw2
 * 
 * The core computes the optimal scores for global and extension alignment of
 * two sequences using dual gap affine cost function and a custom scoring matrix.
 * The core also supports banded alignment and the z-drop heuristics:
 * https://doi.org/10.1093/bioinformatics/bty191
 * 
 * Before running the mothods of this library, you must ensure that the target FPGA 
 * is properly configured with the ksw2_score core through the fpga_configure function.
 * 
 * @file
 */
#ifndef HUGENOMIC_KSW2_SCORE_HPP
#define HUGENOMIC_KSW2_SCORE_HPP

#include <cinttypes>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace hugenomic
{
	/** 
	 * The maximum supported sequence size
	 */
	const uint64_t KSW2S_MAX_SEQUENCE_SIZE = 65536;

    /** 
	 * The maximum supported size for the custom scoring matrix
	 */ 
    const uint64_t KSW2S_MAX_SCORING_MATRIX_SIZE = 64;

	/** 
	 * A special value assigned to a max score returned by the aligner when the 
	 * max score is not found (such as when the alignment is z-dropped)
	 */
	const int32_t KSW2S_MAX_NOT_FOUND = -0x40000000;

	/**
	 * A structure that defines the result of the execution of a ksw2 score alignment on two
	 * pair of sequences
	 */
	typedef struct {

		/** 
		 * Max score of extension alignment:
		 * - query from 0 to [X]
		 * - target from 0 to [Y]
		 * 
		 * can be equal to KSW2S_MAX_NOT_FOUND if not found
		 */
        int32_t extension_max_score;

		/** 
		 * Query end position ([X]) for extension alignment 
		 */
		int32_t extension_query_end_position;

		/** 
		 * Target end position ([Y]) for extension alignment
		 */
        int32_t extension_target_end_position;
        
        /** 
		 * Max score of global query alignment:
		 * - query from 0 to query_len - 1 
		 * - target from 0 to [Y]
		 * 
		 * can be equal to KSW2S_MAX_NOT_FOUND if not found
		 */
        int32_t global_query_max_score;

		/** 
		 * Target end position ([Y]) for global query alignment
		 */		
        int32_t global_query_target_end_position;

		/** 
		 * Max score of global target alignment:
		 * - query from 0 to [X]
		 * - target from 0 to target_len - 1
		 * 
		 * can be equal to KSW2S_MAX_NOT_FOUND if not found
		 */	
        int32_t global_target_max_score;

		/** 
		 * Query end position ([X]) for global target alignment 
		 */ 
        int32_t global_target_query_end_position;

		/** 
		 * Max score of global alignment:
		 * - query from 0 to query_len - 1
		 * - target from 0 to target_len - 1
		 * 
		 * can be equal to KSW2S_MAX_NOT_FOUND if not found
		 */
        int32_t global_max_score;					
	
		/** 
		 * Whether the alignment has been z-dropped
		 */
        bool zdropped;
        
	} ksw2s_response;

	/**
	 * Performs a sequence alignment of one or more pairs of sequences on a target FPGA.
	 * 
	 * The ksw2_score core computes at once the score of 4 different types of alignments.
	 * These alignments are: global alignment, global query alignment, global target alignment, 
	 * extension alignment. 
	 * For each alignment type, each pair of query and target sequences are
	 * aligned from their beginning up to an end position (query, target) that depends on the 
	 * alignment type as follows: 
	 *
	 * - **global alignment** (query_length - 1, target_length - 1)
	 * - **global query alignment** (query_length - 1, ``FREE``)
	 * - **global target alignment** (``FREE``, target_length - 1)
	 * - **extension alignment** (``FREE``, ``FREE``)
	 *
	 * Where ``FREE`` means that the algorithm selects the position that maximizes the alignment score.
	 *
	 * The cost function for gaps is computed as the minimum of two affine gap functions 
	 * using the following formula, where l is the length of the gap:
	 * min{start_gap_cost + l * extend_gap_cost, start_gap_cost_2 + l * extend_gap_cost_2}
	 * 
	 * The cost function associated to a match/mismatch is provided by a scoring matrix.
	 * The scoring matrix is an M*M matrix encoded as a linear array that determines 
	 * the score to assign when comparing an element q of the query to an element t 
	 * of the target sequence.
	 * In details, the score is assigned as:
	 * 
	 * score = scoring_matrix[t * M + q]
	 * 
	 * The value M must be in the range [2, KSW2S_MAX_SCORING_MATRIX_SIZE].
	 * 
	 * Notice that the values stored in the query and target sequences must be 
	 * between 0 and M-1 otherwise an std::invalid_argument will be thrown.
	 * 
	 * Furthermore each sequence cannot be longer than KSW2S_MAX_SEQUENCE_SIZE.
	 *
	 * **PERFORMANCE TIP:**
	 * you should always try to maximize the number of pairs of sequences in order to obtain 
	 * maximum performance.
	 * 
	 * @param fpga_id The identifier of the target FPGA to use
	 *
	 * @param pairs_of_sequences The pairs of query and target sequences to align where
	 * 			the first element of each pair is a query sequence and
	 * 			the second element of each pair is a target sequence
	 *
	 * @param start_gap_cost The cost for starting a gap 
	 * @param extend_gap_cost The cost for extending a gap 
	 * @param start_gap_cost_2 The cost for starting a gap (second gap function)
	 * @param extend_gap_cost_2 The cost for extending a gap (second gap function)
	 *
	 * @param scoring_matrix Scoring matrix of size M*M encoded as a linear vector
	 * 
	 * @param bandwidth The width of the band for banded alignment, use a negative number for full alignment
	 * @param zdrop The z value for the z-drop heuristic, use a negative number to avoid z-drop
	 * 
	 * @return Vector of ksw2s_response containing the alignment result for all pairs of sequences
	 * 
	 * @throw std::invalid_argument In case of invalid sequences or scoreing matrix
	 * @throw invalid_fpga_configuration If the target FPGA is not configured with the ksw2s_core core
	 * @throw invalid_fpga_identifier If an invalid FPGA identifier is specified
	 */
	std::vector<ksw2s_response> ksw2s_batch(
			uint32_t fpga_id,
			const std::vector< std::pair< std::vector<uint8_t>, std::vector<uint8_t> > > &pairs_of_sequences, 
			int8_t start_gap_cost, int8_t extend_gap_cost, int8_t start_gap_cost_2, int8_t extend_gap_cost_2,
			const std::vector< int8_t > scoring_matrix, int32_t bandwidth, int32_t zdrop);
}

#endif // HUGENOMIC_KSW2_SCORE_HPP
