help([[
Salmon 1.10.1

Transcript-level quantification from RNA-seq data.

]])

whatis("Name: salmon")
whatis("Version: 1.10.1")
whatis("Category: bioinformatics")
whatis("Description: Transcript-level RNA-seq quantification")

local software_root = os.getenv("SOFTWARE_ROOT") 
local root = pathJoin(software_root, "salmon", "1.10.1")

prepend_path("PATH", pathJoin(root, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(root, "lib"))
prepend_path("LIBRARY_PATH", pathJoin(root, "lib"))
prepend_path("CPATH", pathJoin(root, "include"))

setenv("SALMON_HOME", root)
