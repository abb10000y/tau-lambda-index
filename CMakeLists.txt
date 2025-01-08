cmake_minimum_required(VERSION 3.10)

project(tau-lambda-index)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# 設定優化選項（與 -O3 一致）
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w -DNDEBUG")

add_executable(tau-lambda-index src/main.cc)
add_executable(gen_mf src/gen_mf.cc)

# 啟用測試
enable_testing()

# 設定 Google Test 的路徑
set(GTEST_INCLUDE_DIR "/opt/homebrew/opt/googletest/include")
set(GTEST_LIBRARY "/opt/homebrew/opt/googletest/lib")

# 加入 Google Test
find_package(GTest REQUIRED)
#find_package(GTest REQUIRED PATHS "/opt/homebrew/opt/googletest")
include_directories(${GTEST_INCLUDE_DIR})
include_directories("/usr/local/include")
include_directories("/usr/local/include/sdsl")

link_directories("/usr/local/lib")
link_directories(${GTEST_LIBRARY})

# 添加測試執行檔案
add_executable(MyTests 
    tests/test_main.cc
    include/util/utility.cc
)

target_include_directories(MyTests PRIVATE
    ${CMAKE_SOURCE_DIR}/include      # 包含符號表和工具的頭文件
    ${CMAKE_SOURCE_DIR}/include/symbol_table  # 如果有必要，也可以指定具體目錄
    ${CMAKE_SOURCE_DIR}/include/util      # 如果需要其他目錄
)

# 連結 Google Test 庫
target_link_libraries(MyTests ${GTEST_LIBRARIES} pthread sdsl)

# 定義測試
add_test(NAME RunTests COMMAND MyTests)