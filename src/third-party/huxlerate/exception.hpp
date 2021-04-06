/** 
 * Copyright (c) 2019, Huxelerate S.r.l.
 * All rights reserved.
 * 
 * @file
 */
#ifndef HUGENOMIC_EXCEPTION_HPP
#define HUGENOMIC_EXCEPTION_HPP

#include <stdexcept>
#include <string>

namespace hugenomic 
{
    /** 
     * Thrown when an invalid FPGA identifier is specified.
     */
    class invalid_fpga_identifier : public std::invalid_argument
    {
        public:
        // constructors
        invalid_fpga_identifier(std::string const & s) : std::invalid_argument{s} {}
        invalid_fpga_identifier() : std::invalid_argument{"Invalid FPGA identifier"} {}
        invalid_fpga_identifier(int fpga_id) : std::invalid_argument{
            std::string("Invalid FPGA identifier. Unable to find FPGA with id: ") + std::to_string(fpga_id)
        } {}
    };

    /**
     * Thrown when a method of a HUGenomic core is called on a FPGA that has not been configured 
     * for that core.
     */
    class invalid_fpga_configuration : public std::invalid_argument
    {
        public:
        // constructors
        invalid_fpga_configuration(std::string const & s) : std::invalid_argument{s} {}
    };

    /** 
     * Thrown when an invalid worker identifier is specified for a HUGenomic core.
     */
    class invalid_worker_identifier : public std::invalid_argument
    {
        public:
        // constructors
        invalid_worker_identifier(std::string const & s) : std::invalid_argument{s} {}
        invalid_worker_identifier() : std::invalid_argument{"Invalid worker identifier"} {}
    };
}

#endif // HUGENOMIC_EXCEPTION_HPP