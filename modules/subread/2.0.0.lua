help([[
Subread 2.0.0

Includes featureCounts for gene/exon-level read counting.

]])

whatis("Name: subread")
whatis("Version: 2.0.0")
whatis("Category: bioinformatics")
whatis("Description: Subread package including featureCounts")

local software_root = os.getenv("SOFTWARE_ROOT")
local root = pathJoin(software_root, "subread", "2.0.0")

prepend_path("PATH", pathJoin(root, "bin"))

setenv("SUBREAD_HOME", root)
setenv("FEATURECOUNTS_HOME", root)
