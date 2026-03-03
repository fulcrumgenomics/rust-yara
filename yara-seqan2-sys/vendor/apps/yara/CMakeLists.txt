# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/yara/CMakeLists.txt
#
# CMakeLists.txt file for yara.
# ===========================================================================

cmake_minimum_required (VERSION 3.12)
project (seqan_apps_yara CXX)
message (STATUS "Configuring apps/yara")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
if (NOT "${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    find_package (OpenMP COMPONENTS CXX)
    find_package (ZLIB)
    find_package (BZip2)
    find_package (SeqAn CONFIG REQUIRED)
endif ()

# Warn if OpenMP was not found.
if (NOT OPENMP_FOUND)
  message (STATUS "  Could NOT find OpenMP 3.0: building without multi-threading")
endif (NOT OPENMP_FOUND)

# Warn if Zlib was not found.
if (NOT ZLIB_FOUND)
  message (STATUS "  Could NOT find zlib: building without gzipped input and BAM output support")
endif (NOT ZLIB_FOUND)

# Warn if BZip2 was not found.
if (NOT BZIP2_FOUND)
  message (STATUS "  Could NOT find bzip2: building without bzipped input support")
endif (NOT BZIP2_FOUND)

# ----------------------------------------------------------------------------
# App-Level Configuration
# ----------------------------------------------------------------------------

set (SEQAN_APP_VERSION "1.0.7")

option (YARA_LARGE_CONTIGS "Set to OFF to disable support for more than 32k contigs or contigs longer than 4Gbp." ON)
if (YARA_LARGE_CONTIGS)
  add_definitions (-DYARA_LARGE_CONTIGS=1)
endif (YARA_LARGE_CONTIGS)

# ----------------------------------------------------------------------------
# GBWT Integration (Optional)
# ----------------------------------------------------------------------------

option(YARA_WITH_GBWT "Enable GBWT-based HLA optimization" OFF)

if (YARA_WITH_GBWT)
    # Require C++17 for GBWT
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    # Require OpenMP
    if (NOT OPENMP_FOUND)
        message(FATAL_ERROR "GBWT requires OpenMP support")
    endif()

    include(ExternalProject)

    # SDSL-lite (vgteam fork)
    set(SDSL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external)
    ExternalProject_Add(
        sdsl_external
        GIT_REPOSITORY https://github.com/vgteam/sdsl-lite.git
        GIT_TAG master
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${SDSL_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_CXX_STANDARD=17
            -DCMAKE_CXX_FLAGS=-fPIC
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DCMAKE_POLICY_VERSION_MINIMUM=3.5
        BUILD_BYPRODUCTS
            ${SDSL_INSTALL_DIR}/lib/libsdsl.a
            ${SDSL_INSTALL_DIR}/lib/libdivsufsort.a
            ${SDSL_INSTALL_DIR}/lib/libdivsufsort64.a
    )

    # GBWT library
    # GBWT's Makefile expects ../sdsl-lite/Make.helper to exist
    # We create a symlink to satisfy this requirement
    set(SDSL_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/sdsl_external-prefix/src/sdsl_external)
    set(GBWT_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gbwt_external-prefix/src/gbwt_external)
    ExternalProject_Add(
        gbwt_external
        GIT_REPOSITORY https://github.com/jltsiren/gbwt.git
        GIT_TAG master
        DEPENDS sdsl_external
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${SDSL_SOURCE_DIR}
            <SOURCE_DIR>/../sdsl-lite
        BUILD_COMMAND ${CMAKE_COMMAND} -E env
            SDSL_DIR=${SDSL_INSTALL_DIR}
            make -C <SOURCE_DIR> -j4
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory
            <SOURCE_DIR>/include ${SDSL_INSTALL_DIR}/include
            COMMAND ${CMAKE_COMMAND} -E copy
            <SOURCE_DIR>/lib/libgbwt.a ${SDSL_INSTALL_DIR}/lib/libgbwt.a
        BUILD_BYPRODUCTS ${SDSL_INSTALL_DIR}/lib/libgbwt.a
        BUILD_IN_SOURCE TRUE
    )

    # Create imported targets
    add_library(SDSL::sdsl STATIC IMPORTED GLOBAL)
    set_target_properties(SDSL::sdsl PROPERTIES
        IMPORTED_LOCATION ${SDSL_INSTALL_DIR}/lib/libsdsl.a
        INTERFACE_INCLUDE_DIRECTORIES ${SDSL_INSTALL_DIR}/include
        INTERFACE_LINK_LIBRARIES
            "${SDSL_INSTALL_DIR}/lib/libdivsufsort.a;${SDSL_INSTALL_DIR}/lib/libdivsufsort64.a"
    )
    add_dependencies(SDSL::sdsl sdsl_external)

    add_library(GBWT::gbwt STATIC IMPORTED GLOBAL)
    set_target_properties(GBWT::gbwt PROPERTIES
        IMPORTED_LOCATION ${SDSL_INSTALL_DIR}/lib/libgbwt.a
        INTERFACE_INCLUDE_DIRECTORIES ${SDSL_INSTALL_DIR}/include
        INTERFACE_LINK_LIBRARIES SDSL::sdsl
    )
    add_dependencies(GBWT::gbwt gbwt_external)

    # GSSW library (Graph Smith-Waterman)
    set(GSSW_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gssw_external-prefix/src/gssw_external)
    ExternalProject_Add(
        gssw_external
        GIT_REPOSITORY https://github.com/vgteam/gssw.git
        GIT_TAG master
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make -C <SOURCE_DIR> -j4
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
            <SOURCE_DIR>/src/gssw.h ${SDSL_INSTALL_DIR}/include/gssw.h
            COMMAND ${CMAKE_COMMAND} -E copy
            <SOURCE_DIR>/lib/libgssw.a ${SDSL_INSTALL_DIR}/lib/libgssw.a
            COMMAND ${CMAKE_COMMAND} -E copy_directory
            <SOURCE_DIR>/src/simde ${SDSL_INSTALL_DIR}/include/simde
        BUILD_BYPRODUCTS ${SDSL_INSTALL_DIR}/lib/libgssw.a
        BUILD_IN_SOURCE TRUE
    )

    add_library(GSSW::gssw STATIC IMPORTED GLOBAL)
    set_target_properties(GSSW::gssw PROPERTIES
        IMPORTED_LOCATION ${SDSL_INSTALL_DIR}/lib/libgssw.a
        INTERFACE_INCLUDE_DIRECTORIES ${SDSL_INSTALL_DIR}/include
    )
    add_dependencies(GSSW::gssw gssw_external)

    add_definitions(-DYARA_WITH_GBWT=1)
    set(GBWT_LIBRARIES GBWT::gbwt GSSW::gssw)
endif()

# ----------------------------------------------------------------------------
# Disable Compiler Warnings
# ----------------------------------------------------------------------------

# Disables warning about conversion from 64 bit to 32 bit integrals.
# Note that there is a runtime exception thrown if the type is not sufficient.
# Unfortunately the VC does not detect it.
if (MSVC)
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4244")  # When compiling on 32 bit windows.
    set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} /wd4267")  # When compiling on 64 bit windows.
endif (MSVC)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add definitions set by the build system.
add_definitions (-DSEQAN_REVISION="${SEQAN_REVISION}")
add_definitions (-DSEQAN_DATE="${SEQAN_DATE}")
add_definitions (-DSEQAN_APP_VERSION="${SEQAN_APP_VERSION}")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Update the list of file names below if you add source files to your application.
add_executable(yara_indexer indexer.cpp
                            misc_tags.h
                            misc_timer.h
                            misc_options.h
                            misc_types.h
                            store_seqs.h
                            bits_matches.h
                            index_fm.h)

add_executable(yara_mapper    mapper.cpp
                              mapper.h
                              mapper_collector.h
                              mapper_classifier.h
                              mapper_ranker.h
                              mapper_filter.h
                              mapper_extender.h
                              mapper_verifier.h
                              mapper_aligner.h
                              mapper_writer.h
                              basic_alphabet.h
                              file_pair.h
                              file_prefetched.h
                              store_seqs.h
                              misc_tags.h
                              misc_timer.h
                              misc_options.h
                              misc_types.h
                              bits_reads.h
                              bits_hits.h
                              bits_matches.h
                              bits_context.h
                              bits_seeds.h
                              bits_bucket.h
                              find_extender.h
                              find_verifier.h
                              index_fm.h)

# Add GBWT headers and sources if enabled
if (YARA_WITH_GBWT)
    target_sources(yara_indexer PRIVATE
                   gbwt_hla_index.h
                   gssw_aligner.h
                   gssw_aligner_impl.h
                   gssw_aligner.cpp)
    target_sources(yara_mapper PRIVATE
                   gbwt_hla_index.h
                   gbwt_hla_finder.h
                   gbwt_hla_extender.h
                   gssw_aligner.h
                   gssw_aligner_impl.h
                   gssw_aligner.cpp)
endif()

# Add dependencies found by find_package (SeqAn).
target_link_libraries (yara_indexer ${SEQAN_LIBRARIES})
target_link_libraries (yara_mapper ${SEQAN_LIBRARIES})

# Add GBWT libraries if enabled
if (YARA_WITH_GBWT)
    target_link_libraries(yara_indexer ${GBWT_LIBRARIES})
    target_link_libraries(yara_mapper ${GBWT_LIBRARIES})
endif()

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (yara)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install yara in ${PREFIX}/bin directory
install (TARGETS yara_indexer yara_mapper
         DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/yara for SeqAn release builds.
install (FILES LICENSE README.rst
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

seqan_add_app_test (yara)

# ----------------------------------------------------------------------------
# Setup Common Tool Description for Generic Workflow Nodes
# ----------------------------------------------------------------------------

# Include executable yara_mapper in CTD structure.
set (SEQAN_CTD_EXECUTABLES ${SEQAN_CTD_EXECUTABLES}
                           yara_indexer yara_mapper CACHE INTERNAL "")

#Add test to check correct ctd generation.
seqan_add_ctd_test (yara_indexer)
seqan_add_ctd_test (yara_mapper)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:yara")
  set (CPACK_PACKAGE_NAME "yara")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Yara - Yet Another Read Aligner")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Enrico Siragusa <enrico.siragusa@fu-berlin.de>")
  set (CPACK_PACKAGE_VENDOR "SeqAn Team, FU Berlin")
  set (CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.rst")

  seqan_configure_cpack_app (yara "yara")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:yara")
