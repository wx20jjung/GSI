help([[
]])

prepend_path("MODULEPATH", "/contrib/spack-stack/spack-stack-1.9.0/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

local python_ver=os.getenv("python_ver") or "3.11.7"
local stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
local stack_impi_ver=os.getenv("stack_mpi_ver") or "2021.13"
local cmake_ver=os.getenv("cmake_ver") or "3.30.2"
local prod_util_ver=os.getenv("prod_util_ver") or "2.1.1"
local spack_apps_ver=os.getenv("spack_apps_ver") or "2024.11"

load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_mpi_ver))
load(pathJoin("python", python_ver))
load(pathJoin("cmake", cmake_ver))
load(pathJoin("spack-apps", spack_apps_ver))

load("gsi_common")
load(pathJoin("prod_util", prod_util_ver))

pushenv("CFLAGS", "-xHOST")
pushenv("FFLAGS", "-xHOST")
pushenv("USE_BUFR4", "YES")

pushenv("GSI_BINARY_SOURCE_DIR", "/scratch1/NCEPDEV/global/glopara/fix/gsi/20241022")

whatis("Description: GSI environment on Hera with Intel Compilers")
