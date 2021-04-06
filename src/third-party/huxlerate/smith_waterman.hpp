/**
 * Copyright (c) 2019, Huxelerate S.r.l.
 * All rights reserved.
 *
 * Provides a set of APIs to interact with the HUGenomic smith_waterman core
 * that performs accelerated sequence alignment on FPGA.
 *
 * The core supports local, global and semi-global (a.k.a. overlap) pairwise alignments
 * with custom affine scoring system.
 * For more information about pairwise alignments refer to the documentation of the Seqan library documentation:
 * https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/Alignment/PairwiseSequenceAlignment.html 
 * 
 * Before running the mothods of this library, you must ensure that the target FPGA is properly configured
 * with the smith_waterman core through the fpga_configure function.
 * 
 * @file
 */
#ifndef HUGENOMIC_SMITH_WATERMAN_HPP
#define HUGENOMIC_SMITH_WATERMAN_HPP

#include <cinttypes>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace hugenomic
{
	/** 
	 * The maximum supported size of the smaller among a pair of horizontal and vertical sequences
	 */
	const uint64_t SW_MAX_SIZE_OF_SMALLER_SEQUENCE = 65535;

	/**
	 * The maximum supported size of a sequence
	 */
	const uint64_t SW_MAX_SEQUENCE_SIZE = 715762091;

	/**
	 * A structure that defines the result of the execution of the smith_waterman core on two
	 * pair of sequences
	 */
	typedef struct {
		/**
		 * the optimal alignment score between the two sequences
		 */
		int32_t alignment_score;

		/**
		 * The end position of the vertical sequence in the optimal alignment
		 */
		uint32_t vertical_sequence_end_position;	

		/**
		 * The end position of the horizontal sequence in the optimal alignment
		 */
		uint32_t horizontal_sequence_end_position;

		/**
		 * The start position of the vertical sequence in the optimal alignment
		 */
		uint32_t vertical_sequence_start_position;

		/**
		 * The start position of the horizontal sequence in the optimal alignment
		 */
		uint32_t horizontal_sequence_start_position;

	} sw_response;

	/**
	 * Performs an optimal sequence alignment of one or more pairs of sequences on a target FPGA.
	 * 
	 * Sequences cannot be longer than SW_MAX_SEQUENCE_SIZE and the smaller sequence of a pair 
	 * cannot exceed SW_MAX_SIZE_OF_SMALLER_SEQUENCE.
	 *
	 * **PERFORMANCE TIP:**
	 * you should always try to maximize the number of pairs of sequences in order to obtain 
	 * maximum performance.
	 *
	 * @param fpga_id The identifier of the target FPGA to use
	 *
	 * @param pairs_of_sequences The pairs of horizontal and vertical sequences to align where
	 * 			the first element of each pair is a horizontal sequence and
	 * 			the second element of each pair is a vertical sequence
	 *
	 * @param match_score The score associated to a match
	 * @param mismatch_score The score associated to a mismatch
	 * @param extend_gap_score The score associated to a gap extension
	 * @param start_gap_score The score associated to starting a gap
	 *
	 * @param global_alignment Whether to perform global alignment (true) or local alignment (false)
	 * @param free_top Whether there is no penalty for skipping characters at the
	 *			beginning of the vertical sequence (only for global alignment)
	 * @param free_left Whether there is no penalty for skipping characters at the
	 * 			beginning of the horizontal sequence (only for global alignment)
	 * @param free_right Whether there is no penalty for skipping characters at the
	 * 			end of the horizontal sequence (only for global alignment)
	 * @param free_down Whether there is no penalty for skipping characters at the
	 * 			end of the vertical sequence (only for global alignment)
	 * 
	 * @return Vector of sw_response containing the alignment result for the pairs of sequences
	 * 
	 * @throw std::invalid_argument if a sequence is empty or a pair of sequences exceeds the size constraints
	 * @throw invalid_fpga_configuration if the target FPGA is not configured with the smith_waterman core
	 * @throw invalid_fpga_identifier if an invalid FPGA identifier is specified
	 */
	std::vector<sw_response> sw_batch(
			uint32_t fpga_id,
			const std::vector< std::pair< std::string, std::string > > &pairs_of_sequences, 
			int32_t match_score, int32_t mismatch_score, int32_t extend_gap_score, int32_t start_gap_score,
			bool global_alignment, bool free_top, bool free_left, bool free_right, bool free_down);

	/**
	 * Runs an all-to-all optimal sequence alignment on a given FPGA.
	 *
	 * The method receives as input two vector of std::string containing the pairs of horizontal and vertical sequences
	 * to align against each other and return a vector of sw_response of size num_horizontal_sequences X num_vertical_sequences 
	 * containing the alignment results for all the pairs of sequences.
	 *
	 * Sequences cannot be longer than SW_MAX_SEQUENCE_SIZE and the smaller sequence of a pair of horizontal 
	 * and vertical sequences cannot exceed SW_MAX_SIZE_OF_SMALLER_SEQUENCE.
	 *
	 * **PERFORMANCE TIP:**
	 * Whenever possible try to maximize the number of alignments to perform in a single function call.
	 *
	 * @param fpga_id The identifier of the target FPGA
	 *
	 * @param horizontal_sequences Vector containing all the horizontal sequences
	 * @param vertical_sequences Vector containing all the vertical sequences
	 *
	 * @param match_score The score associated to a match
	 * @param mismatch_score The score associated to a mismatch
	 * @param extend_gap_score The score associated to a gap extension
	 * @param start_gap_score The score associated to starting a gap
	 *
	 * @param global_alignment Whether to perform global alignment (true) or local alignment (false)
	 * @param free_top Whether there is no penalty for skipping characters at the
	 *					beginning of the vertical sequence (only for global alignment)
	 * @param free_left Whether there is no penalty for skipping characters at the
	 * 					beginning of the horizontal sequence (only for global alignment)
	 * @param free_right Whether there is no penalty for skipping characters at the
	 * 					end of the horizontal sequence (only for global alignment)
	 * @param free_down Whether there is no penalty for skipping characters at the
	 * 					end of the vertical sequence (only for global alignment)
	 * 
	 * @return Vector of sw_response containing all the alignment results stored in row-wise format.
	 * 					row 0: from  result[0]                                to  result[1 * horizontal_sequences.size() - 1]
	 * 					row 1: from  result[1 * horizontal_sequences.size()]  to  result[2 * horizontal_sequences.size() - 1]
	 * 					row 2: from  result[2 * horizontal_sequences.size()]  to  result[3 * horizontal_sequences.size() - 1]
	 * 					...
	 *
	 * @throw std::invalid_argument If there is an empty sequence or a pair of 
	 * 			vertical and horizontal sequences exceed size constraints (see above)
	 * @throw invalid_fpga_configuration If the target FPGA is not configured with the smith_waterman core
	 * @throw invalid_fpga_identifier If an invalid FPGA identifier is specified
	 */
	std::vector<sw_response> sw_all_to_all(uint32_t fpga_id,
			const std::vector<std::string> &horizontal_sequences,
			const std::vector<std::string> &vertical_sequences,
			int32_t match_score, int32_t mismatch_score, int32_t extend_gap_score, int32_t start_gap_score,
			bool global_alignment,
			bool free_top, bool free_left, bool free_right, bool free_down);
}

#endif // HUGENOMIC_SMITH_WATERMAN_HPP
