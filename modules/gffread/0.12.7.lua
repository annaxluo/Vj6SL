help([[
gffread 0.12.7

Utility for processing GFF/GTF annotation files and extracting transcript
sequences.
]])

whatis("Name: gffread")
whatis("Version: 0.12.7")
whatis("Category: bioinformatics")
whatis("Description: GFF/GTF annotation processing utility")

local software_root = os.getenv("SOFTWARE_ROOT")
local root = pathJoin(software_root, "gffread", "0.12.7")

prepend_path("PATH", pathJoin(root, "bin"))

setenv("GFFREAD_HOME", root)
