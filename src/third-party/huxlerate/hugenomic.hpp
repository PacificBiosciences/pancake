/** 
 * Copyright (c) 2019, Huxelerate S.r.l.
 * All rights reserved.
 * 
 * @file
 */
#ifndef HUGENOMIC_HUGENOMIC_HPP
#define HUGENOMIC_HUGENOMIC_HPP

#include "exception.hpp"

// Header files containing the methods available for each of the HUGenomic cores. 
#include "cores/smith_waterman.hpp"
#include "cores/dynamic_time_warping.hpp"
#include "cores/ksw2_score.hpp"
#include "cores/ksw2_cigar.hpp"

namespace hugenomic
{
	
	/** 
	 * The HUGenomic cores available within the HUGenomic library.
	 */
	enum hugenomic_core {
		smith_waterman,
		dynamic_time_warping,
		ksw2_score,
		ksw2_cigar
	};

	/** 
	 * The maximum number of FPGAs per AWS F1 instance.
	 */
	const uint32_t MAX_NUM_FPGA_PER_INSTANCE = 8;

	/**
	 * Configure an HUGenomic core onto one of the available FPGAs.
	 * The number of available FPGAs varies according to the AWS F1 instance being selected
	 * (https://aws.amazon.com/it/ec2/instance-types/f1/):
	 * - f1.2xlarge: 1 FPGA (fpga id: 0)
	 * - f1.4xlarge: 2 FPGAs (fpga id from 0 to 1)
	 * - f1.16xlarge 8 FPGAs (fpga id from 0 to 7)
	 * 
	 * @param core The HUGenomic core to configure on the FPGA
	 * @param fpga_id Identifier of the FPGA
	 * 
	 * @throw invalid_fpga_identifier if an invalid FPGA identifier is specified
	 */
	void fpga_configure(hugenomic_core core, uint32_t fpga_id);

	/**
	 * Release the resources allocated to a specific FPGA.
	 *
	 * @param fpga_id The identifier of the FPGA (from 0 to 7)
	 * 
	 * @throw invalid_fpga_identifier If an invalid FPGA identifier is specified
	 */
	void fpga_release(uint32_t fpga_id);

	/**
	 * Return the number of available FPGAs on the current instance.
	 * @note When using the software emulator, this method always returns MAX_NUM_FPGA_PER_INSTANCE.
	 * 
	 * @return The number of available FPGAs
	 */
	uint32_t get_available_fpgas();
}

#endif // HUGENOMIC_HUGENOMIC_HPP
