# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-src"
  "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-build"
  "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-subbuild/matplotlib-populate-prefix"
  "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-subbuild/matplotlib-populate-prefix/tmp"
  "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-subbuild/matplotlib-populate-prefix/src/matplotlib-populate-stamp"
  "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-subbuild/matplotlib-populate-prefix/src"
  "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-subbuild/matplotlib-populate-prefix/src/matplotlib-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-subbuild/matplotlib-populate-prefix/src/matplotlib-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/sheridan7/Workspace/Laglos/_deps/matplotlib-subbuild/matplotlib-populate-prefix/src/matplotlib-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
