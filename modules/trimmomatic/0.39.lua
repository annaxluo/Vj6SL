help([[
Trimmomatic 0.39

Read trimming tool for Illumina paired-end and single-end data.
]])

whatis("Name: Trimmomatic")
whatis("Version: 0.39")
whatis("Category: read trimming")
whatis("Description: Adapter and quality trimming tool")

local software_root = os.getenv("SOFTWARE_ROOT") 
local root = pathJoin(software_root, "Trimmomatic", "0.39")

setenv("TRIMMOMATIC_HOME", root)
setenv("TRIMMOMATIC_JAR", pathJoin(root, "trimmomatic-0.39.jar"))
setenv("TRIMMOMATIC_ADAPTERS", pathJoin(root, "adapters"))
