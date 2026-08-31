help([[
GATK 4.5.0.0

Genome Analysis Toolkit.
]])

whatis("Name: GATK")
whatis("Version: 4.5.0.0")
whatis("Category: variant analysis")
whatis("Description: Genome Analysis Toolkit")

local software_root = os.getenv("SOFTWARE_ROOT")
local root = pathJoin(software_root, "gatk", "4.5.0.0")
local gatk_home = pathJoin(root, "gatk")
local gatk_env = pathJoin(root, "gatk_env")

prepend_path("PATH", gatk_home)

setenv("GATK_HOME", gatk_home)
setenv("GATK_CONDA_ENV", gatk_env)
setenv("PYTHONNOUSERSITE", "1")
