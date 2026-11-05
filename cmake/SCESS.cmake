# SCESS: command line driver for the quantum chemistry methods, built into build/SCESS

add_executable(SCESS
  ${SCESS_SRC_DIR}/SCESS/Application.cpp
  ${SCESS_SRC_DIR}/SCESS/Driver.cpp
  ${SCESS_SRC_DIR}/SCESS/Driver.h
  ${SCESS_SRC_DIR}/SCESS/TestEinsum.cpp
  ${SCESS_SRC_DIR}/SCESS/Tests.cpp
  ${SCESS_SRC_DIR}/SCESS/Tests.h
)
scess_target_defaults(SCESS)

target_link_libraries(SCESS PRIVATE
  QuantumSCESS
  DahliaTools
)
