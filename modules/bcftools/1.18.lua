help([[
bcftools 1.18

Tools for variant calling and manipulating VCF/BCF files.
]])

whatis("Name: bcftools")
whatis("Version: 1.18")
whatis("Category: bioinformatics")
whatis("Description: Tools for VCF/BCF variant file manipulation")

local software_root = os.getenv("SOFTWARE_ROOT")
local root = pathJoin(software_root, "bcftools", "1.18")

prepend_path("PATH", pathJoin(root, "bin"))
prepend_path("MANPATH", pathJoin(root, "share", "man"))

setenv("BCFTOOLS_HOME", root)
