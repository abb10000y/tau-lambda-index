cmake_minimum_required(VERSION 3.10)

project(tau-lambda-index)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(tau-lambda-index src/main.cc)

# 啟用測試
enable_testing()

# 加入 Google Test
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# 添加測試執行檔案
add_executable(MyTests tests/test.cc)

# 連結 Google Test 庫
target_link_libraries(MyTests ${GTEST_LIBRARIES} pthread)

# 定義測試
add_test(NAME RunTests COMMAND MyTests)