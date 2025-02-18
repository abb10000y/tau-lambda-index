cmake_minimum_required(VERSION 3.10)
project(rmqrmm)

# 設定 C++ 標準與最佳化
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -march=native")

# 設定 source files
set(SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/includes/Basic_rmq.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RMQRMM64.cpp
)

# 建立靜態函式庫（對應 `ar rc rmqrmmBP.a Basic_rmq.o RMQRMM64.o`）
add_library(rmqrmmBP STATIC ${SOURCES})

# 設定輸出目標（確保編譯後得到 `rmqrmmBP.a`）
set_target_properties(rmqrmmBP PROPERTIES OUTPUT_NAME "rmqrmmBP")