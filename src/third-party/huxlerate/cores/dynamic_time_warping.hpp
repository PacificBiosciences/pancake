/**
 * Copyright (c) 2019, Huxelerate S.r.l.
 * All rights reserved.
 *
 * Provides a set of APIs to interact with the HUGenomic dynamic_time_warping core
 * that performs accelerated signal alignment on FPGA.
 *
 * The core performs an accelerated Dynamic Time Warping algorithm and supports both euclidean and 
 * squared euclidean distance functions.
 *
 * Before running the mothods of this library, you must ensure that the target FPGA is properly configured
 * with the dynamic_time_warping core through the fpga_configure function.
 * 
 * @file
 */

#ifndef HUGENOMIC_DYNAMIC_TIME_WARPING_HPP
#define HUGENOMIC_DYNAMIC_TIME_WARPING_HPP

#include <cinttypes>
#include <cstddef>
#include <string>
#include <vector>

namespace hugenomic
{
	/** 
	 * The maximum supported size of the smaller among a pair of horizontal and vertical sequences
	 */
	const uint64_t DTW_MAX_SIZE_OF_SMALLER_SEQUENCE = 65535;

	/** 
	 * The maximum supported size of a sequence
	 */
	const uint64_t DTW_MAX_SEQUENCE_SIZE = 178891323;

	/**
	 * A structure that defines the result of the execution of a Dynamic Time Warping on two
	 * pair of signals.
	 */
	typedef struct {
		/**
		 * the optimal distance score between the two signals
		 */
		float alignment_distance;

		/**
		 * The start position of the vertical sequence in the optimal alignment
		 */
		uint32_t vertical_sequence_start_position;

		/**
		 * The end position of the vertical sequence in the optimal alignment
		 */
		uint32_t vertical_sequence_end_position;

		/**
		 * The start position of the horizontal sequence in the optimal alignment
		 */
		uint32_t horizontal_sequence_start_position;

		/**
		 * The end position of the horizontal sequence in the optimal alignment
		 */
		uint32_t horizontal_sequence_end_position;

	} dtw_response;

	/**
	 * An enum to specify the type of distance function to use to calculate the optimal distance score
	 */
	enum dtw_distance_function {
		dtw_euclidean,
		dtw_squared_euclidean
	};

	/**
	 * Performs an optimal sequence alignment using the Dynamic Time Warping algorithm 
	 * of one or more pairs of sequences on a target FPGA.
	 * 
	 * Sequences cannot be longer than DTW_MAX_SEQUENCE_SIZE and the smaller sequence of a pair 
	 * cannot exceed DTW_MAX_SIZE_OF_SMALLER_SEQUENCE.
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
	 * @param distance_function The distance function to use
	 * @param free_top Whether there is no penalty for skipping elements at the
	 *			beginning of the vertical sequence (only for global alignment)
	 * @param free_left Whether there is no penalty for skipping elements at the
	 * 			beginning of the horizontal sequence (only for global alignment)
	 * @param free_right Whether there is no penalty for skipping elements at the
	 * 			end of the horizontal sequence (only for global alignment)
	 * @param free_down Whether there is no penalty for skipping elements at the
	 * 			end of the vertical sequence (only for global alignment)
	 * 
	 * @return Vector of dtw_response containing the alignment result for all pairs of sequences	 
	 * 
	 * @throw std::invalid_argument If a sequence is empty or a pair of sequences exceeds the size constraints
	 * @throw invalid_fpga_configuration If the target FPGA is not configured with the dynamic_time_warping core
	 * @throw invalid_fpga_identifier If an invalid FPGA identifier is specified
	 */
	std::vector<dtw_response> dtw_batch(
			uint32_t fpga_id,
			const std::vector< std::pair< std::vector<float>, std::vector<float> > > &pairs_of_sequences, 
			dtw_distance_function distance_function,
			bool free_top, bool free_left, bool free_right, bool free_down);

	/**
	 * Runs an all-to-all sequence alignment using the Dynamic Time Warping algorithm on a given FPGA.
	 *
	 * The method receives as input two vectors of std::vector<float> containing the pairs of horizontal 
	 * and vertical sequences to align against each other and returns a vector of dtw_response of size:
	 * (number of horizontal sequences) X (number of vertical_sequences)
	 * containing the alignment results for all the pairs of sequences.
	 *
	 * Sequences cannot be longer than DTW_MAX_SEQUENCE_SIZE and the smaller sequence of a pair of horizontal 
	 * and vertical sequences cannot exceed DTW_MAX_SIZE_OF_SMALLER_SEQUENCE.
	 *
	 * **PERFORMANCE TIP:**
	 * Whenever possible try to maximize the number of alignments to perform in a single function call.
	 *
	 * @param fpga_id The identifier of the target FPGA
	 *
	 * @param horizontal_sequences Vector containing all the horizontal sequences
	 * @param vertical_sequences Vector containing all the vertical sequences
	 *
	 * @param distance_function Whether to use euclidean or squared euclidean as distance function
	 * @param free_top Whether there is no penalty for skipping characters at the
	 *					beginning of the vertical sequence
	 * @param free_left Whether there is no penalty for skipping characters at the
	 * 					beginning of the horizontal sequence
	 * @param free_right Whether there is no penalty for skipping characters at the
	 * 					end of the horizontal sequence
	 * @param free_down Whether there is no penalty for skipping characters at the
	 * 					end of the vertical sequence
	 * 
	 * @return Vector of dtw_response containing all the alignment results stored in row-wise format.
	 * 					row 0: from  result[0]                                to  result[1 * horizontal_sequences.size() - 1]
	 * 					row 1: from  result[1 * horizontal_sequences.size()]  to  result[2 * horizontal_sequences.size() - 1]
	 * 					row 2: from  result[2 * horizontal_sequences.size()]  to  result[3 * horizontal_sequences.size() - 1]
	 * 					...
	 *
	 * @throw std::invalid_argument If there are no horizontal or vertical sequences or if there is a pair 
	 * 			of vertical and horizontal sequences which are either empty or exceed size constraints (see above)
	 * @throw invalid_fpga_configuration If the target FPGA is not configured with the dynamic_time_warping core
	 * @throw invalid_fpga_identifier If an invalid FPGA identifier is specified
	 * @throw std::logic_error If there are pending requests enqueued to any of the worker of the selected FPGA
	 */
	std::vector<dtw_response> dtw_all_to_all(uint32_t fpga_id,
			const std::vector<std::vector<float> > &horizontal_sequences,
			const std::vector<std::vector<float> > &vertical_sequences,
			dtw_distance_function distance_function,
			bool free_top, bool free_left, bool free_right, bool free_down);
}

#endif // HUGENOMIC_DYNAMIC_TIME_WARPING_HPP
