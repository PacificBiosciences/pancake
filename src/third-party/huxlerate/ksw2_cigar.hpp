/**
 * Copyright (c) 2021, Huxelerate S.r.l.
 * All rights reserved.
 *
 * Provides a set of APIs to interact with the HUGenomic ksw2_cigar core
 * that performs accelerated sequence alignment on FPGA.
 *
 * The core functionality is based on the ksw2_extd method from the ksw2 library: 
 * https://github.com/lh3/ksw2
 * 
 * The core computes the optimal scores and the CIGAR for global and extension alignment 
 * of two sequences using dual gap affine cost function and a custom scoring matrix.
 * The core also supports banded alignment and the z-drop heuristics:
 * https://doi.org/10.1093/bioinformatics/bty191
 * 
 * @note If you only need the optimal scores but not the CIGAR it is recommended to use
 * the ksw2_score core since it provides better performance.
 *
 * Before running the mothods of this library, you must ensure that the target FPGA 
 * is properly configured with the ksw2_cigar core through the fpga_configure function.
 * 
 * @file
 */
#ifndef HUGENOMIC_KSW2_CIGAR_HPP
#define HUGENOMIC_KSW2_CIGAR_HPP

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
	const uint64_t KSW2C_MAX_SEQUENCE_SIZE = 65536;

    /** 
	 * The maximum supported size for the custom scoring matrix
	 */ 
    const uint64_t KSW2C_MAX_SCORING_MATRIX_SIZE = 64;

	/** 
	 * A special value assigned to a max score returned by the aligner when the 
	 * max score is not found (such as when the alignment is z-dropped)
	 */
	const int32_t KSW2C_MAX_NOT_FOUND = -0x40000000;

    /**
     * The data type used for storing a CIGAR operation within a CIGAR (e.g. 10M, 2I, 4D).
     * A CIGAR operation encodes both the operation type and its length.
     * @note We use the same encoding of the ksw2 library.
     */
    typedef uint32_t ksw2c_cigar_op;

    /**
     * Retrieves the operation type of a CIGAR operation.
     *
     * @param cigar_op The cigar operation
     * @return The operation type, can be one of:
     *  - 'M': match/mismatch
     *  - 'I': insertion 
     *  - 'D': deletion
     *  - 'U': undefined (only if the cigar operation is incorrect)
     */
    inline char ksw2c_get_cigar_op_type(ksw2c_cigar_op cigar_op) {
        return "MIDU"[cigar_op & 0x3];
    }

    /**
     * Retrieves the length of a CIGAR operation.
     *
     * @param cigar_op The CIGAR operation
     * @return The length of the CIGAR operation 
     */
    inline uint32_t ksw2c_get_cigar_op_length(ksw2c_cigar_op cigar_op) {
        return cigar_op >> 4;
    }

	/**
	 * A structure that defines the result of the execution of a ksw2_cigar alignment on two
	 * pair of sequences
	 */
	typedef struct {

		/** 
		 * Max score of extension alignment:
		 * - query from 0 to [X]
		 * - target from 0 to [Y]
		 * 
		 * can be equal to KSW2C_MAX_NOT_FOUND if not found
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
		 * can be equal to KSW2C_MAX_NOT_FOUND if not found
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
		 * can be equal to KSW2C_MAX_NOT_FOUND if not found
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
		 * can be equal to KSW2C_MAX_NOT_FOUND if not found
		 */
        int32_t global_max_score;	

        /**
         * The CIGAR of one of the following alignments:
         * - global alignment: 
         * if extension_only_cigar=false
         * - global query alignment: 
         * if extension_only_cigar=true and global_query_max_score + end_bonus > extension_max_score
         * - extension alignment: 
         * if extension_only_cigar=true and global_query_max_score + end_bonus <= extension_max_score
         *
         * @note If the alignment has been z-dropped, the CIGAR of extension alignment is returned by default
         */
        std::vector<ksw2c_cigar_op> cigar;

        /**
         * Whether the CIGAR for global query alignment is returned
         */
        bool global_query_cigar;
	
		/** 
		 * Whether the alignment has been z-dropped
		 */
        bool zdropped;
        
	} ksw2c_response;

	/**
	 * Performs a sequence alignment of one or more pairs of query and target sequences 
	 * on a target FPGA and returns the alignment scores and the CIGAR for each pair of sequences.
	 *
	 * The ksw2_cigar core computes at once the score of 4 different types of alignments.
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
	 * The core returns only the CIGAR of one alignment among: extension, global and global query alignment.
	 * To specify which CIGAR to return the user can use the extension_only_cigar and the end_bonus
	 * parameters. The following rules determine which CIGAR is returned:
	 * - **global alignment CIGAR** \n 
	 *   when extension_only_cigar = false
	 * - **global query alignment CIGAR** \n
	 *   when extension_only_cigar = true AND global_query_max_score + end_bonus > extension_max_score
	 * - **extension alignment CIGAR** \n
	 *   when extension_only_cigar = true AND global_query_max_score + end_bonus <= extension_max_score
	 * 
	 * However, if the alignment is z-dropped, the **extension alignment CIGAR** is returned by default.
	 *
	 * The cost function for gaps is computed as the minimum of two affine gap functions 
	 * using the following formula, where l is the length of the gap:
	 *
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
	 * The value M must be in the range [2, KSW2C_MAX_SCORING_MATRIX_SIZE].
	 * 
	 * Notice that the values stored in the query and target sequences must be 
	 * between 0 and M-1 otherwise an std::invalid_argument will be thrown.
	 * 
	 * Furthermore each sequence cannot be longer than KSW2C_MAX_SEQUENCE_SIZE.
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
     * @param end_bonus The bonus score to determine if the CIGAR of global query alignment
	 * 			should be returned instead of the CIGAR for extension alignment
	 *			when extension_only_cigar is true
     * @param extension_only_cigar Whether to return the CIGAR for extension alignment/global query alignment 
	 *			(true) or the CIGAR for global alignment (false). If the alignment is z-dropped, the flag is
	 *			not considered and the CIGAR for extension alignment is returned by default.
	 * 
	 * @return Vector of ksw2c_response containing the alignment result for all pairs of sequences
	 * 
	 * @throw std::invalid_argument In case of invalid sequences or scoring matrix
	 * @throw invalid_fpga_configuration If the target FPGA is not configured with the ksw2c_core
	 * @throw invalid_fpga_identifier If an invalid FPGA identifier is specified
	 */
	std::vector<ksw2c_response> ksw2c_batch(
			uint32_t fpga_id,
			const std::vector< std::pair< std::vector<uint8_t>, std::vector<uint8_t> > > &pairs_of_sequences, 
			int8_t start_gap_cost, int8_t extend_gap_cost, int8_t start_gap_cost_2, int8_t extend_gap_cost_2,
			const std::vector< int8_t > scoring_matrix, int32_t bandwidth, int32_t zdrop, int32_t end_bonus, 
            bool extension_only_cigar);
}

#endif // HUGENOMIC_KSW2_CIGAR_HPP
