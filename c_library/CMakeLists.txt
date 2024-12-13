project(ImagineModels VERSION 0.2 DESCRIPTION "A library for Galactic models.")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")

set(IM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/headers/)
INCLUDE_DIRECTORIES ( "/usr/local/include/eigen3 " )
set(IM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/source)


set(INCLUDE_INSTALL_DIR
"${CMAKE_INSTALL_INCLUDEDIR}/ImagineModels")
# CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where ImagineModel header files are installed"


set(IM_SRC_FILES
  #  ${IM_SOURCE_DIR}/helpers.cc
    ${IM_SOURCE_DIR}/helix.cc
    ${IM_SOURCE_DIR}/jaffe.cc
    ${IM_SOURCE_DIR}/regularjf12.cc
    ${IM_SOURCE_DIR}/ymw16.cc
    ${IM_SOURCE_DIR}/tf17.cc
    ${IM_SOURCE_DIR}/sun.cc
    ${IM_SOURCE_DIR}/pshirkov.cc
    ${IM_SOURCE_DIR}/archimedes.cc
    ${IM_SOURCE_DIR}/han.cc
    ${IM_SOURCE_DIR}/wmap.cc
    ${IM_SOURCE_DIR}/stanevbss.cc
    ${IM_SOURCE_DIR}/fauvet.cc
    ${IM_SOURCE_DIR}/tinyakovtkachev.cc
    ${IM_SOURCE_DIR}/hararimollerachroulet.cc  
    ${IM_SOURCE_DIR}/svt22.cc   
)

if(FFTW_FOUND)
    set(IM_SRC_FILES
        ${IM_SRC_FILES}
        ${IM_SOURCE_DIR}/randomscalarfield.cc
        ${IM_SOURCE_DIR}/randomvectorfield.cc
        ${IM_SOURCE_DIR}/gaussianscalar.cc
        ${IM_SOURCE_DIR}/lognormal.cc
        ${IM_SOURCE_DIR}/randomjf12.cc
        ${IM_SOURCE_DIR}/ensslinsteininger.cc
    )
endif()

add_library(ImagineModels SHARED ${IM_SRC_FILES})
target_include_directories(ImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModels)

if (FFTW_FOUND)
    target_include_directories(ImagineModels PUBLIC ${IM_INCLUDE_DIR}/ImagineModelsRandom)
endif()
target_link_libraries(ImagineModels PUBLIC ${LIBRARIES})
if (autodiff_FOUND) 
    target_include_directories(ImagineModels PUBLIC /usr/local/include/autodiff)
    target_include_directories(ImagineModels PUBLIC /usr/include/eigen3)
    target_link_libraries(ImagineModels PUBLIC autodiff::autodiff)
endif()

set_target_properties(ImagineModels PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION}
)

enable_testing()
set(TESTSOURCES grid parameter_update positions)
foreach(test ${TESTSOURCES})
    set(test_name ${test}_ctest)
    add_executable(${test_name} "${PROJECT_SOURCE_DIR}/test/regular/test_${test}.cc")
    target_link_libraries(${test_name} ImagineModels ${LIBRARIES})
    add_test(${test_name} ${test_name})
endforeach()

configure_file(${PROJECT_SOURCE_DIR}/ImagineModels.pc.in ${CMAKE_BINARY_DIR}/ImagineModels.pc @ONLY)

install(TARGETS ImagineModels LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${CMAKE_BINARY_DIR}/ImagineModels.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

