help([[
]])

prepend_path("MODULEPATH", "/apps/contrib/spack-stack/spack-stack-1.9.0/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

local stack_python_ver=os.getenv("stack_python_ver") or "3.10.8"
local stack_intel_ver=os.getenv("stack_intel_ver") or "2024.2.1"
local stack_intel_mpi_ver=os.getenv("stack_intel_mpi_ver") or "2021.13"
local cmake_ver=os.getenv("cmake_ver") or "3.27.9"
local prod_util_ver=os.getenv("prod_util_ver") or "2.1.1"

load(pathJoin("stack-oneapi", stack_intel_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_intel_mpi_ver))
load(pathJoin("python", stack_python_ver))
load(pathJoin("cmake", cmake_ver))

load("gsi_common")
load(pathJoin("prod_util", prod_util_ver))

pushenv("CFLAGS", "-xHOST")
pushenv("FFLAGS", "-xHOST")
pushenv("USE_BUFR4", "YES")

pushenv("GSI_BINARY_SOURCE_DIR", "/work/noaa/global/glopara/fix/gsi/20241022")

whatis("Description: GSI environment on Orion with Intel Compilers")
