help([[
htslib 1.18

Provides utilities such as bgzip and tabix.

]])

whatis("Name: htslib")
whatis("Version: 1.18")
whatis("Category: bioinformatics")
whatis("Description: HTSlib library and utilities for high-throughput sequencing data")

local software_root = os.getenv("SOFTWARE_ROOT") 
local root = pathJoin(software_root, "htslib", "1.18")

prepend_path("PATH", pathJoin(root, "bin"))
prepend_path("LD_LIBRARY_PATH", pathJoin(root, "lib"))
prepend_path("LIBRARY_PATH", pathJoin(root, "lib"))
prepend_path("CPATH", pathJoin(root, "include"))
prepend_path("PKG_CONFIG_PATH", pathJoin(root, "lib", "pkgconfig"))
prepend_path("MANPATH", pathJoin(root, "share", "man"))

setenv("HTSLIB_HOME", root)
