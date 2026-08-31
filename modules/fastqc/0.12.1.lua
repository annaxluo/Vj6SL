help([[
FastQC 0.12.1

Quality control tool for high-throughput sequencing reads.

]])

whatis("Name: FastQC")
whatis("Version: 0.12.1")
whatis("Category: bioinformatics")
whatis("Description: Quality control for high-throughput sequencing reads")

local software_root = os.getenv("SOFTWARE_ROOT")
local root = pathJoin(software_root, "FastQC", "0.12.1")

prepend_path("PATH", root)

setenv("FASTQC_HOME", root)
