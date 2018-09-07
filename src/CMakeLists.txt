# Parts of code copied from Seqan docs on CMake, and others from kallisto's own
# CMake files.
file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)
list(REMOVE_ITEM sources debugUtil.cpp)

add_library(thing_core ${sources} ${headers})

add_executable(thing main.cpp)
add_executable(debug debugUtil.cpp)

find_package(ZLIB)
find_package(BZip2)
find_package(SeqAn REQUIRED)
include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})

find_package(Threads)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
target_link_libraries(thing thing_core ${SEQAN_LIBRARIES} pthread)
target_link_libraries(debug thing_core)
