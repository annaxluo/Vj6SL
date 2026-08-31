help([[
STAR 2.7.8a

Align Reads to a Reference.
]])

whatis("Name: STAR")
whatis("Version: 2.7.8a")
whatis("Category: aligner")
whatis("Description: Spliced RNA-seq aligner")

local software_root = os.getenv("SOFTWARE_ROOT") 
local root = pathJoin(software_root, "STAR", "2.7.8a")

prepend_path("PATH", pathJoin(root, "bin"))

setenv("STAR_HOME", root)
