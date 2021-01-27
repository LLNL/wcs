#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <string>
#include "wcs_params.hpp"

namespace po = boost::program_options;

namespace wcs {


bool cmd_line_opts::parse_cmd_line(int argc, char** argv)
{
  try {
    po::positional_options_description pdesc;
    pdesc.add("input-model", 1);

    std::string usage
      = std::string("Usage: ") + argv[0] + " [options] input-model\n"
      + "Allowed options";
    po::options_description desc(usage);
    desc.add_options()
      ("help", "Show the help message.")
      ("setup-all", po::value<std::string>(),
        "Specify the merged prototext file for all setups.")
      ("setup-sim", po::value<std::string>(),
        "Specify the prototext file for simulation setup.")
      ("setup-part", po::value<std::string>(),
        "Specify the prototext file for partitioning setup.")
      ("setup-des", po::value<std::string>(),
        "Specify the prototext file for discrete event processing setup.")
      ("input-model", po::value<std::string>(),
        "Specify an input model file. (required)")
      ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(desc).positional(pdesc).run(), vm);

    po::notify(vm);

    m_is_set = false;

    if (vm.size() == 0u || vm.count("help")) {
      std::cout << desc << std::endl;
      return true;
    }

    if (vm.count("input-model")) {
      m_input_model = vm["input-model"].as<std::string>();
      std::cout << "Input model file: "
                << m_input_model << std::endl;
    } else {
      std::cerr << "An input model file is required." << std::endl;
      return false;
    }

    if (vm.count("setup-all")) {
      m_all_setup = vm["setup-all"].as<std::string>();
      m_is_set = true;
      std::cout << "Merged setup file: "
                << m_all_setup << std::endl;
      return true;
    }

    if (vm.count("setup-sim")) {
      m_sim_setup = vm["setup-sim"].as<std::string>();
      m_is_set = true;
      std::cout << "Simulation setup file: "
                << m_sim_setup << std::endl;
    }
    if (vm.count("setup-part")) {
      m_part_setup = vm["setup-part"].as<std::string>();
      m_is_set = true;
      std::cout << "Partition setup file: "
                << m_part_setup << std::endl;
    }
    if (vm.count("setup-des")) {
      m_des_setup = vm["setup-des"].as<std::string>();
      m_is_set = true;
      std::cout << "Discrete event processing setup file: "
                << m_des_setup << std::endl;
    }
    if (!m_all_setup.empty()) {
      m_sim_setup.clear();
      m_part_setup.clear();
      m_des_setup.clear();
      m_is_set = true;
    }
  } catch(std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
      return false;
  } catch(...) {
    std::cerr << "Unknown exception!" << std::endl;
  }

  return true;
}

void cmd_line_opts::show() const
{
  std::string msg;
  msg = "Command line options used:";
  msg += "\n - input model: " + m_input_model;

  if (m_all_setup.empty()) {
    msg += "\n - simulation setup: " + m_sim_setup;
    msg += "\n - partition setup: " + m_part_setup;
    msg += "\n - DES setup: " + m_des_setup;
  } else {
    msg += "\n - all setup: " + m_all_setup;
  }

  std::cout << msg << std::endl << std::endl;
}

} // end of namespace wcs

#if 0 // For testing
int main(int argc, char** argv)
{
  wcs::cmd_line_opts cmd;
  bool ok = cmd.parse_cmd_line(argc, argv);
  if (!ok) return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
#endif
