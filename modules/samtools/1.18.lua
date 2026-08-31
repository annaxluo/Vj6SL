help([[
Samtools 1.18

Tools for manipulating SAM/BAM/CRAM files.
]])

whatis("Name: samtools")
whatis("Version: 1.18")
whatis("Category: bioinformatics")
whatis("Description: Tools for manipulating high-throughput sequencing alignments")

local software_root = os.getenv("SOFTWARE_ROOT") 
local root = pathJoin(software_root, "samtools", "1.18")

prepend_path("PATH", pathJoin(root, "bin"))
prepend_path("MANPATH", pathJoin(root, "share", "man"))

setenv("SAMTOOLS_HOME", root)
