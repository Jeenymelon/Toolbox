add_executable(TestMultiThreadPrinter MultiThreadPrinter_t.cpp)
target_link_libraries(TestMultiThreadPrinter
        -llog4cplus)