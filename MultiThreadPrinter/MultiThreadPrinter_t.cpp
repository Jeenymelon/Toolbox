#include "MultiThreadPrinter.hpp"

void test_monitor() {
  REGISTER("Chinese");
  REGISTER("English");
  REGISTER("Japanese");
  REGISTER("赛马娘");

  INFO("Chinese", "你好");
  DEBUG("English", "hello");
  WARN("Japanese", "こんにちは");
  ERROR("赛马娘", "曼波");
}

int main () {
    test_monitor();
    return 0;
}